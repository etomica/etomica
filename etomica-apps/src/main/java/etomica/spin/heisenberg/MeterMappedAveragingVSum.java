package etomica.spin.heisenberg;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManager;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

import static etomica.math.SpecialFunctions.besselI;

public class MeterMappedAveragingVSum implements IDataSource, AgentSource<MoleculeAgent> {

    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final Space space;
    protected final PotentialCompute potentialMaster;
    protected final NeighborIterator nbrIterator;
    protected PotentialCallbackPhiSumHeisenberg secondDerivativeSum;
    protected PotentialCalculationMoleculeAgentSum vSum;
    protected PotentialCalculationMoleculeAgentSumMinusIdeal vSumMinusIdeal;
    protected PotentialCalculationMoleculeAgentSumPair vSumPair;
    protected PotentialCalculationMoleculeAgentSumMinusIdealPair vSumPairMinusIdeal;
    protected final double temperature;
    protected final double bJ;
    protected final double mu;
    protected final double bt;
    double I0bJ, I1bJ, I2bJ;
    double[] InbJArray, Inm1bJArray, Inm2bJArray, Inp1bJArray;
    protected boolean doIdeal, doPair, doVSum, doVSumMI, doAEEMF;
    protected int nMax;
    protected Vector dr;
    protected Vector tmp, workVector;
    protected final AtomLeafAgentManager<MoleculeAgent> leafAgentManager;
    private final Box box;

    protected final int N;
    protected final MeterMeanField meterMeanField;
    protected double[] nbrSsum = new double[0];//nbrSsum[i] holds sum sin(thetaj) for nbrs j of atom i
    protected double[] nbrCsum = new double[0];//same but cos(thetaj)
    protected final Vector[] dtdotkdt0k, dtdotkdetak;
    protected MeterEnergyMeanField.PotentialCalculationCSsum pcCSsum;

    protected PotentialCalculationHeisenberg Ans;
    protected PotentialCalculationPhiijMF pcPhiIJ;
    protected double phiIJsum;

    public MeterMappedAveragingVSum(Box box, double temperature, double interactionS, double dipoleMagnitude, PotentialCompute potentialMaster, NeighborManager nbrManager, boolean doIdeal, boolean doPair, boolean doVSum, boolean doVSumMI, boolean doAEEMF, int nMax) {
        int nValues = 28;
        data = new DataDoubleArray(nValues);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{nValues});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.box = box;
        this.space = box.getSpace();
        this.temperature = temperature;
        this.potentialMaster = potentialMaster;
        this.nbrIterator = nbrManager.makeNeighborIterator();
        bt = 1 / temperature;
        bJ = interactionS * bt;
        mu = dipoleMagnitude;
        this.doIdeal = doIdeal;
        this.doPair = doPair;
        this.doVSum = doVSum;
        this.doVSumMI = doVSumMI;
        this.doAEEMF = doAEEMF;
        this.nMax = nMax;

        I0bJ = besselI(0, bJ);
        I1bJ = besselI(1, bJ);
        I2bJ = besselI(2, bJ);
        InbJArray = new double[nMax + 1];
        Inm1bJArray = new double[nMax + 1];
        Inm2bJArray = new double[nMax + 1];
        Inp1bJArray = new double[nMax + 1];
        for (int n = 1; n <= nMax; n++) {
            InbJArray[n] = besselI(n, bJ);
            Inm1bJArray[n] = besselI(n - 1, bJ);
            Inm2bJArray[n] = besselI(n - 2, bJ);
            Inp1bJArray[n] = besselI(n + 1, bJ);
        }

        dr = space.makeVector();
        workVector = space.makeVector();
        tmp = space.makeVector();
        leafAgentManager = new AtomLeafAgentManager<>(this, box);
//        FSum = new PotentialCalculationFSum(space, dipoleMagnitude, interactionS, bt);
        secondDerivativeSum = new PotentialCallbackPhiSumHeisenberg(space, leafAgentManager);

        if (doPair)
            Ans = new PotentialCalculationHeisenberg(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);
        if (doVSum)
            vSum = new PotentialCalculationMoleculeAgentSum(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);
        if (doVSumMI)
            vSumMinusIdeal = new PotentialCalculationMoleculeAgentSumMinusIdeal(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);
        if (doVSum)
            vSumPair = new PotentialCalculationMoleculeAgentSumPair(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);
        if (doVSumMI)
            vSumPairMinusIdeal = new PotentialCalculationMoleculeAgentSumMinusIdealPair(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);
        if (doAEEMF)
            pcPhiIJ = new PotentialCalculationPhiijMF();

        N = box.getLeafList().size();
        meterMeanField = doAEEMF ? new MeterMeanField(space, box, interactionS, nbrManager, temperature) : null;
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
        if (box == null) throw new IllegalStateException("no box");
        IAtomList leafList = box.getLeafList();
        secondDerivativeSum.zeroSum();
        potentialMaster.computeAll(true, secondDerivativeSum);
        Vector[] torques = potentialMaster.getTorques();
        for (int i=0; i<leafList.size(); i++) {
            leafAgentManager.getAgent(leafList.get(i)).torque.E(torques[i]);
        }
        if (doPair) {
            Ans.zeroSum();
        }
        if (doVSum) {
            vSum.zeroSum();
        }
        if (doVSum) {
            vSumPair.zeroSum();
        }
        for (IAtom a1 : box.getLeafList()) {
            nbrIterator.iterUpNeighbors(a1.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
                @Override
                public void accept(IAtom jAtom, Vector rij) {
                    if (doPair) {
                        Ans.go((IAtomOriented) a1, (IAtomOriented) jAtom);
                    }
                    if (doVSum) {
                        vSum.go((IAtomOriented) a1, (IAtomOriented) jAtom);
                        vSumPair.go((IAtomOriented) a1, (IAtomOriented) jAtom);
                    }
                }
            });

        }
        double bt2 = bt * bt;
        double mu2 = mu * mu;
        double bmu = bt * mu;
        int nM = leafList.size();


        double AEE = 0, JEMUEx = 0, JEMUEy = 0, JEMUExSquare = 0, JEMUEySquare = 0;
        double JEEMJEJESelf = 0, UEESelf = 0;
        if (doVSum) {
            for (int i = 0; i < nM; i++) {
                MoleculeAgent agentAtomI = leafAgentManager.getAgent(leafList.get(i));

                //-dvEEi/dti
                double dvEEi = agentAtomI.dvEEx().getX(0) + agentAtomI.dvEEy().getX(0);
                JEEMJEJESelf += dvEEi;

                //-vEi*d2vE/dtidti
                double vExi = agentAtomI.vEx().getX(0);
                double vEyi = agentAtomI.vEy().getX(0);
                double d2vExi = agentAtomI.d2vEx().getX(0);
                double d2vEyi = agentAtomI.d2vEy().getX(0);
                JEEMJEJESelf += vExi * d2vExi + vEyi * d2vEyi;

                //-vEEi*fi
                double fi = bt * agentAtomI.torque.getX(0);
                double vEEi = agentAtomI.vEEx().getX(0) + agentAtomI.vEEy().getX(0);
                UEESelf -= vEEi * fi;
                //-fi*dvEi/dti*vEi
                double dvExi = agentAtomI.dvEx().getX(0);
                double dvEyi = agentAtomI.dvEy().getX(0);
                UEESelf -= fi * (dvExi * vExi + dvEyi * vEyi);
                //vEi*phiii*vEi
                double phiii = bt * agentAtomI.phi.component(0, 0);
                UEESelf += vExi * phiii * vExi + vEyi * phiii * vEyi;
                //-2*vEi*fEi
                //fExi = -bmu Sin[thetai]
                //fEyi = bmu Cos[thetai]
                IAtomOriented atom = (IAtomOriented) leafList.get(i);
                double fExi = -bmu * atom.getOrientation().getDirection().getX(1);
                double fEyi = bmu * atom.getOrientation().getDirection().getX(0);
                UEESelf -= 2 * (vExi * fExi + vEyi * fEyi);
                //-var[JEUME]
                double tmp = dvExi + bmu * atom.getOrientation().getDirection().getX(0) + vExi * fi;
                JEMUExSquare += tmp * tmp;
                tmp = dvEyi + bmu * atom.getOrientation().getDirection().getX(1) + vEyi * fi;
                JEMUEySquare += tmp * tmp;

                JEMUEx += dvExi + bmu * atom.getOrientation().getDirection().getX(0) + vExi * fi;
                JEMUEy += dvEyi + bmu * atom.getOrientation().getDirection().getX(1) + vEyi * fi;
            }
            AEE = -vSumPair.getSumJEEMJEJE() - JEEMJEJESelf + UEESelf + vSumPair.getSumUEE() - JEMUEx * JEMUEx - JEMUEy * JEMUEy;
        }

        double torqueScalar = 0;
        if (doIdeal) {
            dr.E(0);
            tmp.E(0);
            for (int i = 0; i < nM; i++) {
                MoleculeAgent agentAtomI = leafAgentManager.getAgent(leafList.get(i));
                torqueScalar = agentAtomI.torque.getX(0);
                IAtomOriented atom = (IAtomOriented) leafList.get(i);
                dr.PEa1Tv1(torqueScalar, atom.getOrientation().getDirection());
                tmp.PE(atom.getOrientation().getDirection());
            }//i loop
        }
        if (doPair) {
            x[0] = -nM * bt2 * mu2 - bt2 * bt2 * mu2 * dr.squared() + bt * bt2 * mu2 * secondDerivativeSum.getSum()
                    - Ans.getSumJEEMJEJE() + Ans.getSumUEE() - Ans.getSumJEMUExSquare() - Ans.getSumJEMUEySquare()
                    - Ans.getAEEJ0();
        }
        if (doIdeal) {
            x[1] = -nM * bt2 * mu2 - bt2 * bt2 * mu2 * dr.squared() + bt * bt2 * mu2 * secondDerivativeSum.getSum();
            x[2] = mu * tmp.getX(0);
            x[3] = mu * tmp.getX(1);
        }

        if (doPair) {
            x[4] = -Ans.getSumJEEMJEJE() + Ans.getSumUEE() - Ans.getSumJEMUExSquare() - Ans.getSumJEMUEySquare();
        }

        if (doVSum) {
            x[5] = AEE;
            x[6] = JEMUEx;
            x[7] = JEMUEy;
            x[8] = JEEMJEJESelf + vSumPair.getSumJEEMJEJE();
            x[9] = UEESelf + vSumPair.getSumUEE();
            x[10] = JEMUEx * JEMUEx;
            x[11] = JEMUEy * JEMUEy;
        }

        if (doVSumMI) {
            for (int i = 0; i < nM; i++) {
                getAns(leafList.get(i));
            }

            tmp.E(0);
            vSumMinusIdeal.zeroSum();
            for (IAtom a1 : box.getLeafList()) {
                nbrIterator.iterUpNeighbors(a1.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
                    @Override
                    public void accept(IAtom jAtom, Vector rij) {
                        vSumMinusIdeal.go((IAtomOriented) a1, (IAtomOriented) jAtom);
                    }
                });
            }

            double bmu2 = bmu * bmu;
            for (int i = 0; i < nM; i++) {
                MoleculeAgent agentAtomI = leafAgentManager.getAgent(leafList.get(i));
                IAtomOriented atom = (IAtomOriented) leafList.get(i);

                dr.E(atom.getOrientation().getDirection());
                tmp.PE(atom.getOrientation().getDirection());

                double sinti = dr.getX(1);
                double costi = dr.getX(0);
                double cos2ti = 2 * costi * costi - 1;

                double vExiIdeal = -bmu * sinti;
                double vEyiIdeal = bmu * costi;
                double vEExiIdeal = 0.5 * bmu2 * costi * sinti;
                double vEEyiIdeal = -0.5 * bmu2 * costi * sinti;
                double dvExidtiIdeal = -bmu * costi;
                double dvEyidtiIdeal = -bmu * sinti;
                double dvEExidtiIdeal = 0.5 * bmu2 * cos2ti;
                double dvEEyidtiIdeal = -dvEExidtiIdeal;
                double d2vExidtidtiIdeal = bmu * sinti;
                double d2vEyidtidtiIdeal = -bmu * costi;

                agentAtomI.vEx().PE(-3 * vExiIdeal);
                agentAtomI.vEy().PE(-3 * vEyiIdeal);
                agentAtomI.vEEx().PE(-3 * vEExiIdeal);
                agentAtomI.vEEy().PE(-3 * vEEyiIdeal);
                agentAtomI.dvEx().PE(-3 * dvExidtiIdeal);
                agentAtomI.dvEy().PE(-3 * dvEyidtiIdeal);
                agentAtomI.dvEEx().PE(-3 * dvEExidtiIdeal);
                agentAtomI.dvEEy().PE(-3 * dvEEyidtiIdeal);
                agentAtomI.d2vEx().PE(-3 * d2vExidtidtiIdeal);
                agentAtomI.d2vEy().PE(-3 * d2vEyidtidtiIdeal);

            }
            AEE = 0;
            JEMUEx = 0;
            JEMUEy = 0;
            JEEMJEJESelf = 0;
            UEESelf = 0;
            vSumPairMinusIdeal.zeroSum();
            for (IAtom a1 : box.getLeafList()) {
                nbrIterator.iterUpNeighbors(a1.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
                    @Override
                    public void accept(IAtom jAtom, Vector rij) {
                        vSumPairMinusIdeal.go((IAtomOriented) a1, (IAtomOriented) jAtom);
                    }
                });
            }
            for (int i = 0; i < nM; i++) {
                MoleculeAgent agentAtomI = leafAgentManager.getAgent(leafList.get(i));

                //-dvEEi/dti
                double dvEEi = agentAtomI.dvEEx().getX(0) + agentAtomI.dvEEy().getX(0);
                JEEMJEJESelf += dvEEi;

                //-vEi*d2vE/dtidti
                double vExi = agentAtomI.vEx().getX(0);
                double vEyi = agentAtomI.vEy().getX(0);
                double d2vExi = agentAtomI.d2vEx().getX(0);
                double d2vEyi = agentAtomI.d2vEy().getX(0);
                JEEMJEJESelf += vExi * d2vExi + vEyi * d2vEyi;


                //-vEEi*fi
                double fi = bt * agentAtomI.torque.getX(0);
                double vEEi = agentAtomI.vEEx().getX(0) + agentAtomI.vEEy().getX(0);
                UEESelf -= vEEi * fi;
                //-fi*dvEi/dti*vEi
                double dvExi = agentAtomI.dvEx().getX(0);
                double dvEyi = agentAtomI.dvEy().getX(0);
                UEESelf -= fi * (dvExi * vExi + dvEyi * vEyi);
                //vEi*phiii*vEi
                double phiii = bt * agentAtomI.phi.component(0, 0);
                UEESelf += vExi * phiii * vExi + vEyi * phiii * vEyi;
                //-2*vEi*fEi
                IAtomOriented atom = (IAtomOriented) leafList.get(i);
                double fExi = -bmu * atom.getOrientation().getDirection().getX(1);
                double fEyi = bmu * atom.getOrientation().getDirection().getX(0);
                UEESelf -= 2 * (vExi * fExi + vEyi * fEyi);
                //-var[JEUME]
                JEMUEx += dvExi + bmu * atom.getOrientation().getDirection().getX(0) + vExi * fi;
                JEMUEy += dvEyi + bmu * atom.getOrientation().getDirection().getX(1) + vEyi * fi;

            }
            AEE = -vSumPairMinusIdeal.getSumJEEMJEJE() - JEEMJEJESelf + UEESelf + vSumPairMinusIdeal.getSumUEE() - JEMUEx * JEMUEx - JEMUEy * JEMUEy;


            x[12] = AEE;
            x[13] = JEMUEx;
            x[14] = JEMUEy;
            x[15] = JEEMJEJESelf + vSumPairMinusIdeal.getSumJEEMJEJE();
            x[16] = UEESelf + vSumPairMinusIdeal.getSumUEE();
            x[17] = JEMUEx * JEMUEx;
            x[18] = JEMUEy * JEMUEy;
            x[19] = mu * tmp.getX(0);
            x[20] = mu * tmp.getX(1);


            x[21] = vSumPairMinusIdeal.getSumU_Map();
            x[22] = x[21] * x[21];

            x[23] = vSumPairMinusIdeal.getSumU_Con();
            x[24] = x[23] * x[23];
        }

        if(doAEEMF) {
            meterMeanField.getData();
            pcCSsum.reset();
            for (IAtom a1 : box.getLeafList()) {
                nbrIterator.iterUpNeighbors(a1.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
                    @Override
                    public void accept(IAtom jAtom, Vector rij) {
                        pcCSsum.go((IAtomOriented) a1, (IAtomOriented) jAtom);
                    }
                });
            }

            MeterMappedAveragingCorrelation.computeTdotDerivs(bt, meterMeanField, dtdotkdetak, dtdotkdt0k);

            x[25] = 0.0;
            x[26] = 0.0;
            for(int k  = 0; k < N; k++) {
                Vector s = meterMeanField.getSpin(box.getLeafList().get(k));
                x[25] += s.getX(0);
                x[26] += s.getX(1);
            }
            x[27] = x[25] * x[25] + x[26] * x[26]; //variance contribution

            workVector.E(0.0);
            for (IAtom a1 : box.getLeafList()) {
                nbrIterator.iterUpNeighbors(a1.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
                    @Override
                    public void accept(IAtom jAtom, Vector rij) {
                        pcPhiIJ.go((IAtomOriented) a1, (IAtomOriented) jAtom);
                    }
                });
            }

        }

        return data;
    }


    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public MoleculeAgent makeAgent(IAtom a, Box box) {
        return new MoleculeAgent(nMax);
    }

    public void releaseAgent(MoleculeAgent agent, IAtom a, Box box) {

    }

    public void getAns(IAtom atom) {
        MoleculeAgent agentAtom = leafAgentManager.getAgent(atom);
        double bmu = bt * mu;
        double bmu2 = bmu * bmu;

        Vector ei = space.makeVector();
        ei.E(((IAtomOriented) atom).getOrientation().getDirection());
        double t1 = Math.atan2(ei.getX(1), ei.getX(0));


        double cost1 = ei.getX(0);
        double sint1 = ei.getX(1);
        double cos2t1 = 2 * cost1 * cost1 - 1;
        double sin2t1 = 2 * sint1 * cost1;


        for (int n = 0; n <= nMax; n++) {
            if (n == 0) {
                agentAtom.sinntheta[n] = 0;
                agentAtom.cosntheta[n] = 1;
            }

            if (n == 1) {
                agentAtom.sinntheta[n] = sint1;
                agentAtom.cosntheta[n] = cost1;
            }

            if (n == 2) {
                agentAtom.sinntheta[n] = sin2t1;
                agentAtom.cosntheta[n] = cos2t1;
            }

            if (n == 3) {
                agentAtom.sinntheta[n] = sin2t1 * cost1 + cos2t1 * sint1;
                agentAtom.cosntheta[n] = cos2t1 * cost1 - sin2t1 * sint1;
            }

            if (n == 4) {
                agentAtom.sinntheta[n] = sin2t1 * cos2t1 + cos2t1 * sin2t1;
                agentAtom.cosntheta[n] = cos2t1 * cos2t1 - sin2t1 * sin2t1;
            }

            if (n == 5) {
                agentAtom.sinntheta[n] = agentAtom.sinntheta[3] * cos2t1 + agentAtom.cosntheta[3] * sin2t1;
                agentAtom.cosntheta[n] = agentAtom.cosntheta[3] * cos2t1 - agentAtom.sinntheta[3] * sin2t1;
            }


            if (n == 0) {
                agentAtom.Axc0[n] = bmu * (I0bJ + I1bJ) * (-1 + cost1);
                agentAtom.dAxc0[n] = -bmu * (I0bJ + I1bJ) * sint1;
                agentAtom.Axs0[n] = 0;
                agentAtom.dAxs0[n] = 0;
                agentAtom.Axc1[n] = -0.25 * bmu2 * sint1 * sint1 * (I0bJ + 2 * I1bJ + I2bJ);
                agentAtom.dAxc1[n] = -0.25 * bmu2 * (I0bJ + 2 * I1bJ + I2bJ) * sin2t1;
                agentAtom.Axs1[n] = 0;
                agentAtom.dAxs1[n] = 0;

                agentAtom.d2Axc0[n] = -bmu * cost1 * (I0bJ + I1bJ);
                agentAtom.d2Axs0[n] = 0;
                agentAtom.d3Axc0[n] = bmu * sint1 * (I0bJ + I1bJ);
                agentAtom.d3Axs0[n] = 0;
                agentAtom.d2Axc1[n] = -0.5 * bmu2 * cos2t1 * (I0bJ + 2 * I1bJ + I2bJ);
                agentAtom.d2Axs1[n] = 0;


                //y direction
                agentAtom.Ayc0[n] = bmu * (I0bJ + I1bJ) * sint1;
                agentAtom.dAyc0[n] = bmu * (I0bJ + I1bJ) * cost1;
                agentAtom.Ays0[n] = 0;
                agentAtom.dAys0[n] = 0;
                agentAtom.Ayc1[n] = 0.25 * bmu2 * sint1 * sint1 * (I0bJ + 2 * I1bJ + I2bJ);
                agentAtom.dAyc1[n] = 0.25 * bmu2 * sin2t1 * (I0bJ + 2 * I1bJ + I2bJ);
                agentAtom.Ays1[n] = 0;
                agentAtom.dAys1[n] = 0;

                agentAtom.d2Ayc0[n] = -bmu * sint1 * (I0bJ + I1bJ);
                agentAtom.d2Ays0[n] = 0;
                agentAtom.d3Ayc0[n] = -bmu * cost1 * (I0bJ + I1bJ);
                agentAtom.d3Ays0[n] = 0;
                agentAtom.d2Ayc1[n] = 0.5 * bmu2 * (I0bJ + 2 * I1bJ + I2bJ) * cos2t1;
                agentAtom.d2Ays1[n] = 0;


            }


            if (n > 0) {
                int n2 = n * n;
                int n3 = n2 * n;
                int n4 = n2 * n2;
                double InbJ = InbJArray[n];
                double Inm1bJ = Inm1bJArray[n];
                double Inm2bJ = Inm2bJArray[n];
                double Inp1bJ = Inp1bJArray[n];
                double sinnt1 = agentAtom.sinntheta[n];
                double cosnt1 = agentAtom.cosntheta[n];


                double sinnm1t1 = sinnt1 * cost1 - cosnt1 * sint1;
                double sinnp1t1 = sinnt1 * cost1 + cosnt1 * sint1;
                double coshnt1 = Math.cosh(n * t1);
                double sinnm2t1 = sinnt1 * cos2t1 - cosnt1 * sin2t1;
                double sinnp2t1 = sinnt1 * cos2t1 + cosnt1 * sin2t1;
                double cosnm2t1 = cosnt1 * cos2t1 + sinnt1 * sin2t1;
                double cosnp2t1 = cosnt1 * cos2t1 - sinnt1 * sin2t1;
                double cosnm1t1 = cosnt1 * cost1 + sinnt1 * sint1;
                double cosnp1t1 = cosnt1 * cost1 - sinnt1 * sint1;


                agentAtom.Axc0[n] = 2 * bmu * (((bJ + 2 * bJ * n2) * Inm1bJ + (bJ - n + 2 * (1 + bJ) * n2 - 2 * n3) * InbJ) * cost1 * cosnt1
                        + (2 * bJ * Inm1bJ + (1 + 2 * bJ - 2 * n + 2 * n2) * InbJ) * n * sint1 * sinnt1)
                        / (bJ + 4 * bJ * n4);


                agentAtom.dAxc0[n] = (2 * InbJ * (-cosnt1 * sint1 + (n - 2 * n3) * cost1 * sinnt1)
                        + Inm1bJ * (1 + n - 2 * n3) * sinnm1t1
                        + (-1 + n - 2 * n3) * Inp1bJ * sinnp1t1
                ) * bmu / (1 + 4 * n4);

                agentAtom.Axs0[n] = ((1 + 2 * n + 2 * n2) * Inm1bJ * sinnm1t1
                        + 2 * InbJ * (-2 * n * cosnt1 * sint1 + (1 + 2 * n2) * cost1 * sinnt1)
                        + (1 - 2 * n + 2 * n2) * Inp1bJ * sinnp1t1
                ) * bmu / (1 + 4 * n4);

                agentAtom.dAxs0[n] = (-2 * bmu * n * (InbJ + (-1 + 2 * n2) * (-bJ * Inm1bJ + (-bJ + n) * InbJ)) * cost1 * cosnt1
                        - 2 * bmu * (bJ * Inm1bJ + (bJ - n + n2 - 2 * n4) * InbJ) * sint1 * sinnt1
                ) / (bJ + 4 * bJ * n4);

                agentAtom.Axc1[n] = (Inm2bJ * (cosnm2t1 - coshnt1) * n2 / (2 - 2 * n + n2)
                        + (-4 * bJ * (4 + n4) * cosnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                        + 2 * n2 * (2 + n2 - 2 * n) * I0bJ * cosnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                        + bJ * n2 * (2 + 2 * n + n2) * I0bJ * (bJ * InbJ * (cosnm2t1 + coshnt1) + 2 * Inm1bJ * (bJ * cosnm2t1 + (n - 1) * coshnt1))
                ) / (bJ * bJ * (4 + n4) * I0bJ)
                ) * bmu2 / 4.0 / n2;

                agentAtom.dAxc1[n] = 0.25 * bmu2 * (
                        -(n - 2) * Inm2bJ * sinnm2t1 / (2 - 2 * n + n2)
                                +
                                (bJ * (-bJ * n * (-4 - 2 * n + n3) * I0bJ * sinnm2t1 * (2 * Inm1bJ + InbJ) + (16 + 4 * n4) * sinnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                                        - 2 * n * (4 - 2 * n + n3) * I0bJ * sinnp2t1 * (bJ * (bJ - 1 - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)

                                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );

                agentAtom.Axs1[n] = ((n2 * Inm2bJ * sinnm2t1) / (2 - 2 * n + n2)
                        + (bJ * (bJ * n2 * (2 + 2 * n + n2) * I0bJ * sinnm2t1 * (2 * Inm1bJ + InbJ) - (16 + 4 * n4) * sinnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                        + 2 * n2 * (2 - 2 * n + n2) * I0bJ * sinnp2t1 * ((-bJ + bJ * bJ - n * bJ) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n2 + 2 * n) * InbJ)
                ) / (bJ * bJ * (4 + n4) * I0bJ)
                ) * bmu2 / 4.0 / n2;

                agentAtom.dAxs1[n] = 0.25 * bmu2 * (((n - 2) * Inm2bJ * cosnm2t1) / (2 - 2 * n + n2)
                        + (bJ * (bJ * n * (-4 - 2 * n + n3) * I0bJ * cosnm2t1 * (2 * Inm1bJ + InbJ)
                        - 4 * (4 + n4) * cosnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                )
                        + 2 * n * (4 - 2 * n + n3) * I0bJ * cosnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );


                agentAtom.d2Axc0[n] = (-(1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) * Inp1bJ * cosnp1t1
                        + (-1 + n) * Inm1bJ * ((1 + n) * cosnm1t1 - 2 * n3 * cosnm1t1)
                        + InbJ * (-2 * (1 - n2 + 2 * n4) * cost1 * cosnt1 + 4 * n3 * sint1 * sinnt1)
                ) * bmu / (1 + 4 * n4);

                agentAtom.d2Axs0[n] = (-2 * bJ * bmu * Inm1bJ * (2 * n3 * cosnt1 * sint1 + (1 - n2 + 2 * n4) * cost1 * sinnt1)
                        + 2 * bmu * InbJ * (n * (1 - n2 - 2 * bJ * n2 + 2 * n3 + 2 * n4) * cosnt1 * sint1 + (n * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) + bJ * (-1 + n2 - 2 * n4)) * cost1 * sinnt1)
                ) / (bJ + 4 * bJ * n4);

                agentAtom.d3Axc0[n] = (2 * InbJ * ((1 - n2 + 4 * n4) * cosnt1 * sint1 + n * (1 + n2 + 2 * n4) * cost1 * sinnt1)
                        + (1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) * Inp1bJ * sinnp1t1
                        - (-1 + n) * (-1 + n) * Inm1bJ * ((1 + n) * sinnm1t1 + 2 * n3 * sinnm1t1)
                ) * bmu / (1.0 + 4 * n4);

                agentAtom.d3Axs0[n] = -(bJ * Inm1bJ * (n * (1 + n2 + 2 * n4) * cost1 * cosnt1 + (-1 + n2 - 4 * n4) * sint1 * sinnt1)
                        + InbJ * (-n * cost1 * cosnt1 * ((1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) - bJ * (1 + n2 + 2 * n4)) + sint1 * sinnt1 * (n * (1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) + bJ * (-1 + n2 - 4 * n4)))
                ) * 2 * bmu / (bJ + 4 * bJ * n4);

                agentAtom.d2Axc1[n] = 0.25 * bmu2 * (-(n - 2) * (n - 2) * Inm2bJ * cosnm2t1 / (2 - 2 * n + n2)
                        + (bJ * (-bJ * (n - 2) * (n - 2) * n * (2 + 2 * n + n2) * I0bJ * cosnm2t1 * (2 * Inm1bJ + InbJ) + 4 * n * (4 + n4) * cosnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                        - 2 * n * (2 + n) * (4 - 2 * n + n3) * I0bJ * cosnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );

                agentAtom.d2Axs1[n] = 0.25 * bmu2 * (-((n - 2) * (n - 2) * Inm2bJ * sinnm2t1) / (2 - 2 * n + n2)
                        + (bJ * (-bJ * (n - 2) * (n - 2) * n * (2 + 2 * n + n2) * I0bJ * sinnm2t1 * (2 * Inm1bJ + InbJ) + 4 * n * (4 + n4) * sinnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                        - 2 * n * (n + 2) * (4 - 2 * n + n3) * I0bJ * sinnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );

                agentAtom.Ayc0[n] = (2 * bmu * cosnt1 * sint1 * ((bJ + 2 * bJ * n2) * Inm1bJ + (bJ - n + 2 * n2 + 2 * bJ * n2 - 2 * n3) * InbJ)
                        - 2 * bmu * n * cost1 * sinnt1 * (2 * bJ * Inm1bJ + (1 + 2 * bJ + 2 * (-1 + n) * n) * InbJ)
                ) / (bJ + 4 * bJ * n4);

                agentAtom.dAyc0[n] = 2.0 * bmu * ((bJ * Inm1bJ + (bJ - n + n2 - 2 * n4) * InbJ) * cost1 * cosnt1
                        + n * sint1 * sinnt1 * (InbJ + (-1 + 2 * n2) * (-bJ * Inm1bJ + (n - bJ) * InbJ))
                ) / (bJ + 4 * bJ * n4);

                agentAtom.Ays0[n] = (bJ * (1 - 2 * n + 2 * n2) * Inp1bJ * (-cosnp1t1 + coshnt1)
                        + bJ * Inm1bJ * ((1 + 2 * n + 2 * n2) * cosnm1t1 + (-1 + 2 * n - 2 * n2) * coshnt1)
                        + 2 * InbJ * (2 * bJ * n * cost1 * cosnt1 + n * (1 - 2 * n + 2 * n2) * coshnt1 + bJ * (1 + 2 * n2) * sint1 * sinnt1)
                ) * bmu / (bJ + 4 * bJ * n4);

                agentAtom.dAys0[n] = ((1 + n - 2 * n3) * Inm1bJ * sinnm1t1
                        + 2 * InbJ * (n * (-1 + 2 * n2) * cosnt1 * sint1 + cost1 * sinnt1)
                        + (1 - n + 2 * n3) * Inp1bJ * sinnp1t1
                ) * bmu / (1.0 + 4 * n4);


                agentAtom.Ayc1[n] = (n2 * Inm2bJ * (-cosnm2t1 + coshnt1) / (2 - 2 * n + n2)
                        - (4 * bJ * (4 + n4) * cosnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                        + 2 * n2 * (2 - 2 * n + n2) * I0bJ * cosnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                        + bJ * n2 * (2 + 2 * n + n2) * I0bJ * (bJ * InbJ * (cosnm2t1 + coshnt1) + 2 * Inm1bJ * (bJ * cosnm2t1 + (-1 + n) * coshnt1))
                ) / (bJ * bJ * (4 + n4) * I0bJ)
                ) * bmu2 / 4 / n2;

                agentAtom.dAyc1[n] = 0.25 * bmu2 * (((-2 + n) * Inm2bJ * sinnm2t1) / (2 - 2 * n + n2)
                        + (bJ * (bJ * n * (-4 - 2 * n + n3) * I0bJ * sinnm2t1 * (2 * Inm1bJ + InbJ) + (16 + 4 * n4) * sinnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                        + 2 * n * (4 - 2 * n + n3) * I0bJ * sinnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );

                agentAtom.Ays1[n] = (-bJ * n2 * (2 + n * (2 + n)) * I0bJ * sinnm2t1 * ((-1 + bJ + n) * Inm1bJ + bJ * InbJ)
                        - 2 * bJ * (4 + n4) * sinnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                        - n2 * (2 - 2 * n + n2) * I0bJ * sinnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n * (1 + n)) * InbJ)
                ) * bmu2 / (2 * bJ * bJ * n2 * (4 + n4) * I0bJ);

                agentAtom.dAys1[n] = 0.25 * bmu2 * (-(-2 + n) * Inm2bJ * cosnm2t1 / (2 - 2 * n + n2)
                        + (bJ * (-bJ * n * (-4 - 2 * n + n3) * I0bJ * cosnm2t1 * (2 * Inm1bJ + InbJ) - 4 * (4 + n4) * cosnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                        - 2 * n * (4 - 2 * n + n3) * I0bJ * cosnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );


                agentAtom.d2Ayc0[n] = 2 * bmu * (-(bJ * Inm1bJ + (bJ + n * (-1 + n - 2 * n3)) * InbJ) * cosnt1 * sint1 + n2 * (InbJ + (-1 + 2 * n2) * (-bJ * Inm1bJ + (-bJ + n) * InbJ)) * cosnt1 * sint1 - n * (bJ * Inm1bJ + (bJ + n * (-1 + n - 2 * n3)) * InbJ) * cost1 * sinnt1 + n * (InbJ + (-1 + 2 * n2) * (-bJ * Inm1bJ + (-bJ + n) * InbJ)) * cost1 * sinnt1) / (bJ + 4 * bJ * n4);

                agentAtom.d2Ays0[n] = (bmu * ((-1 + n) * (1 + n - 2 * n3) * Inm1bJ * cosnm1t1 + (1 + n) * (1 - n + 2 * n3) * Inp1bJ * cosnp1t1 + 2 * InbJ * (n * cost1 * cosnt1 + n * (-1 + 2 * n2) * cost1 * cosnt1 - sint1 * sinnt1 - n2 * (-1 + 2 * n2) * sint1 * sinnt1))) / (1 + 4 * n4);


                agentAtom.d3Ayc0[n] = 2 * bmu * (bJ * Inm1bJ * (-(1 - n2 + 4 * n4) * cost1 * cosnt1 + n * (1 + n2 + 2 * n4) * sint1 * sinnt1)
                        + InbJ * (-(-n * (1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) + bJ * (1 - n2 + 4 * n4)) * cost1 * cosnt1 + n * (-(1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) + bJ * (1 + n2 + 2 * n4)) * sint1 * sinnt1)
                ) / (bJ + 4 * bJ * n4);


                agentAtom.d3Ays0[n] = -bmu * (-(-1 + n) * (-1 + n) * (-1 + n) * (1 + 2 * n + 2 * n2) * Inm1bJ * sinnm1t1
                        + 2 * InbJ * (n * (1 + n2 + 2 * n4) * cosnt1 * sint1 + (1 - n2 + 4 * n4) * cost1 * sinnt1)
                        + (1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) * Inp1bJ * sinnp1t1
                ) / (1 + 4 * n4);

                agentAtom.d2Ayc1[n] = 0.25 * bmu2 * ((-2 + n) * (-2 + n) * Inm2bJ * cosnm2t1 / (2 - 2 * n + n2)
                        + (bJ * (bJ * (-2 + n) * (-2 + n) * n * (2 + 2 * n + n2) * I0bJ * cosnm2t1 * (2 * Inm1bJ + InbJ)
                        + 4 * n * (4 + n4) * cosnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                )
                        + 2 * n * (2 + n) * (4 - 2 * n + n3) * I0bJ * cosnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );


                agentAtom.d2Ays1[n] = 0.25 * bmu2 * ((n - 2) * (n - 2) * Inm2bJ * sinnm2t1 / (2 - 2 * n + n2)
                        + (bJ * (bJ * (n - 2) * (n - 2) * n * (2 + 2 * n + n2) * I0bJ * sinnm2t1 * (2 * Inm1bJ + InbJ)
                        + 4 * n * (4 + n4) * sinnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                )
                        + 2 * n * (2 + n) * (4 - 2 * n + n3) * I0bJ * sinnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );


            }
        }

    }

    private class PotentialCalculationPhiijMF {
        private final Vector myWorkVector = space.makeVector();

        public void go(IAtomOriented iatom, IAtomOriented jatom) {
            Vector io = iatom.getOrientation().getDirection();
            Vector jo = jatom.getOrientation().getDirection();
            myWorkVector.E(meterMeanField.getThetaDot(iatom));
            myWorkVector.ME(meterMeanField.getThetaDot(jatom));
            myWorkVector.TE(myWorkVector);
            workVector.PEa1Tv1(bJ * io.dot(jo), myWorkVector);//(thetaDot_i - thetaDot_j)^2 phi_ij
        }
    }


}
