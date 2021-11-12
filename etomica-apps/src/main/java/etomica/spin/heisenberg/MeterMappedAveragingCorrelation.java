package etomica.spin.heisenberg;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.integrator.Integrator;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManager;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space1d.Vector1D;
import etomica.units.dimensions.Null;

public class MeterMappedAveragingCorrelation implements IDataSource, AtomLeafAgentManager.AgentSource<MeterMappedAveragingCorrelation.CorrelationAgent> {
    protected final DataDoubleArray data;
    protected final DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected double J, bJ, J2;
    protected double bt;
    protected final int L;
    protected final int N;
    protected final int arraySize;
    protected double temperature, beta;
    private Box box;
    protected final Space space;
    protected final PotentialCompute potentialMaster;
    protected final NeighborIterator nbrIterator;
    protected final int[] offset;
    protected final int formula;
    protected final int[] nPairs;
    protected AtomLeafAgentManager<CorrelationAgent> leafAgentManager;
    protected final MeterMeanField meterMeanField;
    protected MeterEnergyMeanField.PotentialCalculationCSsum pcCSsum;
    protected double[] nbrSsum = new double[0];//nbrSsum[i] holds sum sin(thetaj) for nbrs j of atom i
    protected double[] nbrCsum = new double[0];//same but cos(thetaj)
    protected final Vector[] dtdotkdt0k, dtdotkdetak;


    public MeterMappedAveragingCorrelation(Box box, double temperature, double interactionS, PotentialCompute potentialMaster, NeighborManager nbrManager, int formula) {
        this.box = box;
        this.space = box.getSpace();
        this.temperature = temperature;
        beta = 1 / temperature;
        this.potentialMaster = potentialMaster;
        nbrIterator = nbrManager.makeNeighborIterator();
        this.formula = formula;
        bt = 1 / temperature;
        J = interactionS;
        bJ = bt * J;
        J2 = J * J;
        N = box.getLeafList().size();
        L = (int) Math.round(Math.sqrt(N));
        int distance = L / 2 + 1;
//        nPairs = new int[-1 + (distance + 1) * distance / 2];
//        data = new DataDoubleArray(-1 + (distance + 1) * distance / 2);
//        dataInfo = new DataDoubleArray.DataInfoDoubleArray("CIJ", Null.DIMENSION, new int[]{-1 + (distance + 1) * distance / 2});
        arraySize = -1 + (distance + 1) * distance / 2;
        nPairs = new int[3 * arraySize];
        data = new DataDoubleArray(3 * arraySize);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("CIJ", Null.DIMENSION, new int[]{3 * arraySize});
        tag = new DataTag();
        dataInfo.addTag(tag);
        offset = new int[distance];
        for (int i = 0, o = -1; i < distance; i++) {
            offset[i] = o - i;
            o += (distance - i);
        }

        leafAgentManager = new AtomLeafAgentManager<CorrelationAgent>(this, box);
        meterMeanField = formula == 3 ? new MeterMeanField(space, box, interactionS, nbrManager, temperature) : null;

        nbrCsum = new double[N];
        nbrSsum = new double[N];
        pcCSsum = new MeterEnergyMeanField.PotentialCalculationCSsum(nbrCsum, nbrSsum);

        dtdotkdt0k = new Vector[N];
        dtdotkdetak = new Vector[N];
        for(int k=0; k<N; k++) {
            dtdotkdetak[k] = space.makeVector();
            dtdotkdt0k[k] = space.makeVector();
        }

    }

    public IData getData() {
        double[] x = data.getData();
        int distance = L / 2 + 1;
        potentialMaster.computeAll(true);
        Vector[] torques = potentialMaster.getTorques();
        for (int i=0; i<box.getLeafList().size(); i++) {
            leafAgentManager.getAgent(box.getLeafList().get(i)).torque.E(torques[i]);
        }
        IAtomList leafList = box.getLeafList();
        for (int i = 0; i < x.length; i++) {
            x[i] = 0;
            nPairs[i] = 0;
        }
        if (formula == 3) {
            //set up meterMeanField to provide terms needed for further calculation
            meterMeanField.getData();

            //compute terms needed to evaluate d(eta_k)/d(theta_i) and d(theta0_k)/d(theta_i) for i nbr of k
            pcCSsum.reset();
            for (IAtom a1 : box.getLeafList()) {
                IAtomOriented atom1 = (IAtomOriented) a1;
                nbrIterator.iterUpNeighbors(a1.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
                    @Override
                    public void accept(IAtom jAtom, Vector rij) {
                        pcCSsum.go(atom1, (IAtomOriented) jAtom);
                    }
                });
            }

            computeTdotDerivs(beta, meterMeanField, dtdotkdetak, dtdotkdt0k);

        }

        //loop over pairs j and k
        for (int j = 0; j < N; j++) {
            int jRow = j / L, jCol = j % L;
            IAtom aj = leafList.get(j);
            CorrelationAgent agentAtomJ = leafAgentManager.getAgent(aj);
            double bfj = bt * agentAtomJ.torque().getX(0);
            double costj = ((IAtomOriented) leafList.get(j)).getOrientation().getDirection().getX(0);
            double sintj = ((IAtomOriented) leafList.get(j)).getOrientation().getDirection().getX(1);
            Vector sj = formula == 3 ? meterMeanField.getSpin(aj) : null;
            Vector tDotj = formula == 3 ? meterMeanField.getThetaDot(aj) : null;
            Vector sumk = formula == 3 ? space.makeVector() : null;
            Vector workVector = formula == 3 ? space.makeVector() : null;


            for (int k = 0; k < j; k++) {
                int kRow = k / L, kCol = k % L;
                IAtom ak = leafList.get(k);
                CorrelationAgent agentAtomK = leafAgentManager.getAgent(ak);
                double bfk = bt * agentAtomK.torque().getX(0);
                double costk = ((IAtomOriented) leafList.get(k)).getOrientation().getDirection().getX(0);
                double sintk = ((IAtomOriented) leafList.get(k)).getOrientation().getDirection().getX(1);
                int JMKRow = jRow - kRow;
                int JMKCol = jCol - kCol;
                if (JMKRow >= distance) JMKRow = L - JMKRow;
                if (JMKCol >= distance) {
                    JMKCol = L - JMKCol;
                } else if (JMKCol <= -distance) {
                    JMKCol += L;
                } else if (JMKCol < 0) {
                    JMKCol *= -1;
                }

                if (JMKCol < JMKRow) {
                    int tmp = JMKCol;
                    JMKCol = JMKRow;
                    JMKRow = tmp;
                }
                int index = offset[JMKRow] + JMKCol;
                if (formula == 0) {//conventional
                    x[index] += costj * costk + sintj * sintk;
                } else if (formula == 1) { //effective-field, first derivative, high temperature
                    x[index] -= 0.5 * (bfj - bfk) * (sintj * costk - costj * sintk);
                } else if (formula == 2) {
                    x[index] += bfj * bfk * (sintj * sintk + costj * costk);
                    if (index == 0) {
                        x[index] += bJ * (costj * costk + sintj * sintk) * (costj + sintj) * (costk + sintk);
                    }
                } else { //effective field, second derivative
                    Vector sk = meterMeanField.getSpin(ak);
                    x[index] += sj.dot(sk);//covariance contribution
                    if (index == 0) { //neighboring spins
                        // thetaDotj thetaDotk phi_jk contribution
                        sumk.Ea1Tv1(bJ * (costj * costk + sintj * sintk), meterMeanField.getThetaDot(ak));//-beta*phiij = +beta J cos(tj - tk)

                        double etak = meterMeanField.getEta()[k];
                        double sindthetak = meterMeanField.getSindtheta()[k];
                        double cosdthetak = meterMeanField.getCosdtheta()[k];
                        double I2etc = meterMeanField.I2etc[k];
                        double I1I0 = meterMeanField.I1I0[k];
                        double cost0 = meterMeanField.getCost0()[k];
                        double sint0 = meterMeanField.getSint0()[k];
                        Vector thetaDotk = meterMeanField.getThetaDot(ak);

                        double bdEtakdthetaj = 0.25 * beta * J2 * (-sintj * nbrCsum[k] + costj * nbrSsum[k]) / etak;
                        double dtheta0kdthetaj = 0.25 * J2 * (costj * nbrCsum[k] + sintj * nbrSsum[k])/ (etak * etak);

                        workVector.setX(0,I2etc * cost0);
                        workVector.setX(1,I2etc * sint0);

                        workVector.PEa1Tv1(-etak * sindthetak, dtdotkdetak[k]);
                        workVector.PEa1Tv1(sindthetak, thetaDotk);

                        sumk.PEa1Tv1(bdEtakdthetaj, workVector);

                        workVector.setX(0,-I1I0 * sint0);
                        workVector.setX(1, I1I0 * cost0);

                        workVector.PEa1Tv1(-beta * etak * sindthetak, dtdotkdt0k[k]);
                        workVector.PEa1Tv1(-beta * etak * cosdthetak, thetaDotk);

                        sumk.PEa1Tv1(dtheta0kdthetaj, workVector);

                        x[index] += sumk.dot(tDotj);

                    }
                }
                nPairs[index] += 1;
                x[index + arraySize] = JMKRow;
                x[index + 2 * arraySize] = JMKCol;
            }
        }
//        for (int i = 0; i < x.length; i++) {
        for (int i = 0; i < arraySize; i++) {
            x[i] /= nPairs[i];
        }
        return data;
    }

    //Computes derivatives of thetaDot_k with respect to eta_k and theta0_k for all spins
    //The meterMeanField passed as an argument should be set up by invoking its getData() method before calling this
    //dtdotkdetak and dtdotkdt0k should be instantiated and filled with Vector2D instances prior to calling
    public static void computeTdotDerivs(double beta, MeterMeanField meterMeanField, Vector[] dtdotkdetak, Vector[] dtdotkdt0k) {
        int N = dtdotkdetak.length;
        for(int k=0; k < N; k++) {
            double dtheta = meterMeanField.getDtheta()[k];
            double cosdtheta = meterMeanField.getCosdtheta()[k];
            double sindtheta = meterMeanField.getSindtheta()[k];
            double t0 = meterMeanField.getTheta0()[k];
            double cost0 = meterMeanField.getCost0()[k];
            double sint0 = meterMeanField.getSint0()[k];
            double I1I0 = meterMeanField.getI1I0()[k];
            double I2etc = meterMeanField.getI2etc()[k];
            double bh = beta * meterMeanField.getEta()[k];
            double pInv = Math.exp(-bh * cosdtheta);

            // Cmn is integral of (cos(t-t0))^m (sin(t-t0))^n Exp[bh cos(t-t0)] for t = t0 to theta
            // Cmnc is integral from 0 to theta
            // Cmns is integral from Pi/2 to theta
            // integrals are computed using functions that return integral of cosx^m sinx^n Exp[bh cosx] for x = 0 to [first argument]
            double C00 = FunctionCosIntegral.cosInt(dtheta, bh);
            double C00c = C00 - FunctionCosIntegral.cosInt(0 - t0, bh);
            double C00s = C00 - FunctionCosIntegral.cosInt(Math.PI/2. - t0, bh);

            double C10 = FunctionCosIntegral.coscosInt(dtheta, bh);
            double C10c = C10 - FunctionCosIntegral.coscosInt(0 - t0, bh);
            double C10s = C10 - FunctionCosIntegral.coscosInt(Math.PI/2. - t0, bh);

            double C20 = FunctionCosIntegral.cos2cosInt(dtheta, bh);
            double C20c = C20 - FunctionCosIntegral.cos2cosInt(0 - t0, bh);
            double C20s = C20 - FunctionCosIntegral.cos2cosInt(Math.PI/2. - t0, bh);

            double C11 = FunctionCosIntegral.sincoscosInt(dtheta, bh);
            double C11c = C11 - FunctionCosIntegral.sincoscosInt(0 - t0, bh);
            double C11s = C11 - FunctionCosIntegral.sincoscosInt(Math.PI/2. - t0, bh);

            double C01 = FunctionCosIntegral.sincosInt(dtheta, bh);
            double C01c = C01 - FunctionCosIntegral.sincosInt(0 - t0, bh);
            double C01s = C01 - FunctionCosIntegral.sincosInt(Math.PI/2. - t0, bh);

            double C02c = C00c - C20c;//integral of sin^2 is integral of 1 - cos^2
            double C02s = C00s - C20s;

            //compute d(thetaDot_k)/d(eta_k)
            double c = -I2etc * C00c * cost0;
            double s = -I2etc * C00s * sint0;

            c += cost0 * C20c - sint0 * C11c;
            s += cost0 * C11s + sint0 * C20s;

            c -= cost0 * I1I0 * C10c;
            s -= sint0 * I1I0 * C10s;

            c -= cosdtheta * (cost0 * C10c - sint0 * C01c);
            s -= cosdtheta * (cost0 * C01s + sint0 * C10s);

            c += cosdtheta * I1I0 * cost0 * C00c;
            s += cosdtheta * I1I0 * sint0 * C00s;

            dtdotkdetak[k].setX(0, c);
            dtdotkdetak[k].setX(1, s);
            dtdotkdetak[k].TE(-beta * pInv);

            //compute d(thetaDot_k)/d(theta0_k)
            c = + I1I0 * sint0 * C00c;
            s = - I1I0 * cost0 * C00s;

            c += bh * (cost0 * C11c - sint0 * C02c);
            s += bh * (cost0 * C02s + sint0 * C11s);

            c -= bh * cost0 * I1I0 * C01c;
            s -= bh * sint0 * I1I0 * C01s;

            c -= bh * sindtheta * (cost0 * C10c - sint0 * C01c);
            s -= bh * sindtheta * (cost0 * C01s + sint0 * C10s);

            c += bh * sindtheta * I1I0 * cost0 * C00c;
            s += bh * sindtheta * I1I0 * sint0 * C00s;

            dtdotkdt0k[k].setX(0, c);
            dtdotkdt0k[k].setX(1, s);
            dtdotkdt0k[k].TE(-pInv);
        }
    }

    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public CorrelationAgent makeAgent(IAtom a, Box agentBox) {

        return new CorrelationAgent(space);
    }

    public void releaseAgent(CorrelationAgent agent, IAtom atom, Box agentBox) {

    }

    public static class CorrelationAgent implements Integrator.Torquable, Integrator.Forcible {  //need public so to use with instanceof
        public final Vector torque;
        public final Vector force;

        public CorrelationAgent(Space space) {
            torque = new Vector1D();
            force = space.makeVector();
        }

        public Vector torque() {
            return torque;
        }

        public Vector force() {
            return force;
        }
    }
}
