package etomica.spin.heisenberg;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManager;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;
import etomica.util.numerical.BesselFunction;

//public class MeterEnergyMeanField implements IDataSource, AtomLeafAgentManager.AgentSource<MeterEnergyMeanField.ForceTorque> {
public class MeterEnergyMeanField implements IDataSource {

    protected final double J, J2;
    protected final DataDoubleArray data;
    protected final DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
//    protected final PotentialCalculationTorqueSum torqueSum;
    protected final NeighborIterator nbrIterator;
    protected final Box box;
    protected final double beta;
//    protected final AtomLeafAgentManager<ForceTorque> torqueAgentManager;
    protected final PotentialCalculationPhiij pcExtra;
    protected PotentialCalculationCSsum pcCSsum;
    protected final PotentialCalculationCvij pcCvij;
    protected double[] bThetadot = new double[0];
    protected double cvSumExtra;
    protected double cvij;
    protected double[] nbrSsum = new double[0];//nbrSsum[i] holds sum sin(thetaj) for nbrs j of atom i
    protected double[] nbrCsum = new double[0];//same but cos(thetaj)

    protected double[] dJdeta = new double[0];
    protected double[] bdJdtheta0 = new double[0];
    //    protected double[] dtdotdeta = new double[0];
//    protected double[] dtdotdtheta0 = new double[0];
    protected double[] eta = new double[0];
    protected final Vector[] netOrientation;

    public MeterEnergyMeanField(Space space, Box box, double J, NeighborManager nbrManager, double temperature) {
        nbrIterator = nbrManager.makeNeighborIterator();
        beta = 1 / temperature;
        this.J = J;
        J2 = J * J;

        pcExtra = new PotentialCalculationPhiij();
        pcCvij = new PotentialCalculationCvij();
//        torqueSum = new PotentialCalculationTorqueSum();
//        torqueAgentManager = new AtomLeafAgentManager<ForceTorque>(this, box, ForceTorque.class);
//        torqueSum.setAgentManager(torqueAgentManager);
        data = new DataDoubleArray(2);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{2});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.box = box;
        netOrientation = new Vector[box.getLeafList().size()];
        for (int i=0; i<netOrientation.length; i++) {
            netOrientation[i] = box.getSpace().makeVector();
        }
    }

    @Override
    public IData getData() {
        IAtomList atoms = box.getLeafList();
        int atomCount = atoms.size();
        if (bThetadot.length != atomCount) {
            bThetadot = new double[atomCount];
            nbrCsum = new double[atomCount];
            nbrSsum = new double[atomCount];
            pcCSsum = new PotentialCalculationCSsum(nbrCsum, nbrSsum);
            dJdeta = new double[atomCount];
            bdJdtheta0 = new double[atomCount];
            eta = new double[atomCount];
        }
        for (int i=0; i<netOrientation.length; i++) {
            netOrientation[i].E(0);
        }
        pcCSsum.reset();
        for (IAtom a1 : box.getLeafList()) {
            IAtomOriented atom1 = (IAtomOriented) a1;
            nbrIterator.iterUpNeighbors(a1.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
                @Override
                public void accept(IAtom jAtom, Vector rij, int n) {
                    IAtomOriented atom2 = (IAtomOriented) jAtom;
                    Vector ei = atom1.getOrientation().getDirection();
                    Vector ej = atom2.getOrientation().getDirection();

                    netOrientation[atom1.getLeafIndex()].PEa1Tv1(J/2, ej);
                    netOrientation[atom2.getLeafIndex()].PEa1Tv1(J/2, ei);

                    pcCSsum.go(atom1, atom2);
                }
            });
        }

        data.E(0);
        double sum = 0;
        double cvSum = 0;
        for (int i = 0; i < atoms.size(); i++) {
            IAtomOriented a = (IAtomOriented) atoms.get(i);
            Vector h = netOrientation[i];
            double hmag = Math.sqrt(h.squared());
            eta[i] = hmag;
            double bh = beta * hmag;
            h.TE(1.0 / hmag);

            double I0 = BesselFunction.I(0, bh);
            double I1 = BesselFunction.I(1, bh);
            double I2 = BesselFunction.I(2, bh);
            double I1I0 = I1 / I0;
            double I2etc = 0.5 * (I2 / I0 - 2 * I1I0 * I1I0 + 1);

            Vector o = a.getOrientation().getDirection();
            double sindtheta = h.getX(0) * o.getX(1) - o.getX(0) * h.getX(1);
            double cosdtheta = h.dot(o);
            double dtheta = Math.atan2(sindtheta, cosdtheta);
            double pInv = Math.exp(-bh * cosdtheta);

            double C0 = FunctionCosIntegral.cosInt(dtheta, bh);
            double C1 = FunctionCosIntegral.coscosInt(dtheta, bh);
            double C2 = FunctionCosIntegral.cos2cosInt(dtheta, bh);
            double SC0 = FunctionCosIntegral.sincosInt(dtheta, bh);
            double SC1 = FunctionCosIntegral.sincoscosInt(dtheta, bh);

            // beta * thetaDot^beta
            double bvb = bh * (C0 * I1I0 - C1) * pInv;
            bThetadot[i] = bvb;

//            sum += hmag  * (-x + cosdtheta - v * (f + hmag * sindtheta) / temperature);
//            sum += hmag  * (-x - v * (f + hmag * sindtheta) / temperature);
            double ui = hmag * (-I1I0 + bvb * sindtheta);
            sum += ui;
//            usum+= hmag * cosdtheta;

            // beta^2 * thetaDot^beta_beta
            double bbvbb = bh * bh * ((C1 * (I1I0 + cosdtheta) - C0 * I1I0 * cosdtheta - C2)
                    + I2etc * C0) * pInv;

            double bJbeta = bh * (I1I0 - cosdtheta + bvb * sindtheta);

//            cvSum += bh * bh * 0.5 * (I2 / I0 - 2 * I1I0 * I1I0 + 1);
//            cvSum += -bh * bvbb * sindtheta;
//            cvSum += -2 * bh * bv * sindtheta;
//            cvSum += -bh * bv * bJbeta * sindtheta;
//            cvSum += -bv * bv * bh * cosdtheta;

            double term1 = bh * bh * I2etc;
            double term2 = -bh * bbvbb * sindtheta;
            double term3 = -2 * bh * bvb * sindtheta;
            double term4 = -bh * bvb * bJbeta * sindtheta;
            double term5 = -bvb * bvb * bh * cosdtheta;

            cvSum += term1 + term2 + term3 + term4 + term5;

            //ij terms for cv
            dJdeta[i] = I1I0 - cosdtheta + bh * I2etc + bvb * sindtheta;
            bdJdtheta0[i] = -bh * (sindtheta + bvb * cosdtheta);
            double bhdtdotdeta = (bbvbb + bvb);
            double dtdotdtheta0 = bh * pInv * (-Math.exp(bh) * (I1I0 - 1)
                    + bh * (SC0 * I1I0 - C0 * sindtheta * I1I0 - SC1 + C1 * sindtheta)
                    - SC0);
            dJdeta[i] += -sindtheta * bhdtdotdeta;
            bdJdtheta0[i] += -bh * sindtheta * dtdotdtheta0;

        }

        //handle phi_ij contribution
        cvSumExtra = 0;
        cvij = 0;
        for (IAtom a1 : box.getLeafList()) {
            IAtomOriented atom1 = (IAtomOriented) a1;
            nbrIterator.iterUpNeighbors(a1.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
                @Override
                public void accept(IAtom jAtom, Vector rij, int n) {
                    IAtomOriented atom2 = (IAtomOriented) jAtom;
                    pcExtra.go(atom1, atom2);
                    pcCvij.go(atom1, atom2);
                }
            });
        }
        cvSum += 2 * beta * J * cvSumExtra;// 2 is to count each ij twice for phi_ij multiplication

//        PotentialCalculationNbrCounter pcCount = new PotentialCalculationNbrCounter();
//        potentialMaster.calculate(box, allAtoms, pcCount);
//        System.out.println("Neighbor count: "+pcCount.nbrCount+", compare to: "+4*atomCount/2);

        double[] y = data.getData();
        y[0] = sum;
        y[1] = beta * beta * sum * sum + cvSum + cvij;
        return data;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    @Override
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    private class PotentialCalculationPhiij {
        public void go(IAtomOriented iatom, IAtomOriented jatom) {
            Vector io = iatom.getOrientation().getDirection();
            Vector jo = jatom.getOrientation().getDirection();
            cvSumExtra += io.dot(jo) * bThetadot[iatom.getLeafIndex()] * bThetadot[jatom.getLeafIndex()];
        }
    }

    /**
     * Used to compute and store sums of cos(theta) and sin(theta) over neighbors of each atom.
     * Needed for mapped-averaging calculation of heat capacity.
     */
    public static class PotentialCalculationCSsum {
        protected final double[] nbrCsum, nbrSsum;
        public PotentialCalculationCSsum(double[] nbrCsum, double[] nbrSsum) {
            this.nbrCsum = nbrCsum;
            this.nbrSsum = nbrSsum;
        }

        public void go(IAtomOriented iatom, IAtomOriented jatom) {
            int i = iatom.getLeafIndex();
            int j = jatom.getLeafIndex();
            Vector io = iatom.getOrientation().getDirection();
            Vector jo = jatom.getOrientation().getDirection();
            nbrCsum[i] += jo.getX(0);
            nbrCsum[j] += io.getX(0);
            nbrSsum[i] += jo.getX(1);
            nbrSsum[j] += io.getX(1);
        }

        public void reset() {
            for (int i = 0; i < nbrCsum.length; i++) {
                nbrCsum[i] = 0;
                nbrSsum[i] = 0;
            }
        }
    }

    private class PotentialCalculationCvij {

        public void go(IAtomOriented atom1, IAtomOriented atom2) {
            cvij += ijContribution(atom1, atom2) + ijContribution(atom2, atom1);
        }

        private double ijContribution(IAtomOriented iatom, IAtomOriented jatom) {
            int i = iatom.getLeafIndex();
            int j = jatom.getLeafIndex();
            Vector io = iatom.getOrientation().getDirection();
            double costi = io.getX(0);
            double sinti = io.getX(1);

            double bdEtajdThetai = 0.25 * beta * J2 * (-sinti * nbrCsum[j] + costi * nbrSsum[j]) / eta[j];
            double dtheta0jdthetai = 0.25 * J2 * (costi * nbrCsum[j] + sinti * nbrSsum[j])/ (eta[j] * eta[j]);

            return bThetadot[i] * (bdJdtheta0[j] * dtheta0jdthetai + dJdeta[j] * bdEtajdThetai);
        }
    }

//    public class PotentialCalculationNbrCounter implements PotentialCalculation {
//        int nbrCount = 0;
//
//        public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
//            nbrCount++;
//        }
//    }
}