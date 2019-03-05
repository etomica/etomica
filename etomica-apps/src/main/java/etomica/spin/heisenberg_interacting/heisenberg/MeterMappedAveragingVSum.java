package etomica.spin.heisenberg_interacting.heisenberg;

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
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterMappedAveragingVSum implements IDataSource, AgentSource<MoleculeAgent> {
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final Space space;
    protected final PotentialMaster potentialMaster;
    protected final IteratorDirective allAtoms;
    protected PotentialCalculationEnergySum energySum;
    //    protected PotentialCalculationFSum FSum;
    protected PotentialCalculationTorqueSum torqueSum;
    protected PotentialCalculationPhiSum secondDerivativeSum;
    protected PotentialCalculationPhiSumHeisenberg secondDerivativeSumIdeal;
    protected PotentialCalculationMoleculeAgentSum vSum;
    protected PotentialCalculationMoleculeAgentSumMinusIdeal vSumMinusIdeal;
    protected PotentialCalculationMoleculeAgentSumPair vSumPair;
    protected PotentialCalculationMoleculeAgentSumMinusIdealPair vSumPairMinusIdeal;
//    protected PotentialCalculationPair pair;
    protected double temperature;
    protected double J;
    protected double mu;
    protected double bt;
    protected boolean doIdeal, doPair, doVSum, doVSumMI;
    protected Vector dr, tmp;
    protected Vector work;
    protected AtomLeafAgentManager<MoleculeAgent> leafAgentManager;
    private Box box;
    protected PotentialCalculationHeisenberg Ans;

    public MeterMappedAveragingVSum(final Space space, Box box, Simulation sim, double temperature, double interactionS, double dipoleMagnitude, PotentialMaster potentialMaster, boolean doIdeal, boolean doPair, boolean doVSum, boolean doVSumMI) {
        int nValues = 19;
        data = new DataDoubleArray(nValues);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{nValues});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.box = box;
        this.space = space;
        this.temperature = temperature;
        this.potentialMaster = potentialMaster;
        J = interactionS;
        bt = 1 / temperature;
        mu = dipoleMagnitude;
        this.doIdeal = doIdeal;
        this.doPair = doPair;
        this.doVSum = doVSum;
        this.doVSumMI = doVSumMI;

        dr = space.makeVector();
        work = space.makeVector();
        tmp = space.makeVector();
        leafAgentManager = new AtomLeafAgentManager<MoleculeAgent>(this, box, MoleculeAgent.class);
        torqueSum = new PotentialCalculationTorqueSum();
        torqueSum.setAgentManager(leafAgentManager);
//        FSum = new PotentialCalculationFSum(space, dipoleMagnitude, interactionS, bt);
        energySum = new PotentialCalculationEnergySum();
        secondDerivativeSum = new PotentialCalculationPhiSum();
        secondDerivativeSum.setAgentManager(leafAgentManager);
        int nMax = 10;


        if (doIdeal) secondDerivativeSumIdeal = new PotentialCalculationPhiSumHeisenberg(space);

        if (doPair)
            Ans = new PotentialCalculationHeisenberg(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);
        if (doVSum)
            vSum = new PotentialCalculationMoleculeAgentSum(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);
        if (doVSumMI)
            vSumMinusIdeal = new PotentialCalculationMoleculeAgentSumMinusIdeal(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);
//        pair = new PotentialCalculationPair(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);
        if (doVSum)
            vSumPair = new PotentialCalculationMoleculeAgentSumPair(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);
        if (doVSumMI)
            vSumPairMinusIdeal = new PotentialCalculationMoleculeAgentSumMinusIdealPair(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);
        allAtoms = new IteratorDirective();

    }

    public IData getData() {
        double[] x = data.getData();
        if (box == null) throw new IllegalStateException("no box");
        IAtomList leafList = box.getLeafList();
        torqueSum.reset();
        potentialMaster.calculate(box, allAtoms, torqueSum);

        secondDerivativeSum.reset();
        potentialMaster.calculate(box, allAtoms, secondDerivativeSum);

        if (doPair) {
            Ans.zeroSum();
            potentialMaster.calculate(box, allAtoms, Ans);
        }
        if (doVSum) {
            vSum.zeroSum();
            potentialMaster.calculate(box, allAtoms, vSum);
        }
//        pair.zeroSum();
//        potentialMaster.calculate(box, allAtoms, pair);

        if (doVSum) {
            vSumPair.zeroSum();
            potentialMaster.calculate(box, allAtoms, vSumPair);
        }

        if (doIdeal) {
            secondDerivativeSumIdeal.zeroSum();
            potentialMaster.calculate(box, allAtoms, secondDerivativeSumIdeal);
        }
        double bt2 = bt * bt;
        double mu2 = mu * mu;
        double bmu = bt * mu;
        int nM = leafList.getAtomCount();


        double AEE = 0, JEMUEx = 0, JEMUEy = 0, JEMUExSquare = 0, JEMUEySquare = 0;
        double JEEMJEJESelf = 0, UEESelf = 0;
        if (doVSum) {
            for (int i = 0; i < nM; i++) {
                MoleculeAgent agentAtomI = leafAgentManager.getAgent(leafList.getAtom(i));

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
                IAtomOriented atom = (IAtomOriented) leafList.getAtom(i);
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
                MoleculeAgent agentAtomI = leafAgentManager.getAgent(leafList.getAtom(i));
                torqueScalar = agentAtomI.torque.getX(0);
                IAtomOriented atom = (IAtomOriented) leafList.getAtom(i);
                dr.PEa1Tv1(torqueScalar, atom.getOrientation().getDirection());
                tmp.PE(atom.getOrientation().getDirection());
            }//i loop
        }
        if (doPair) {
            x[0] = -nM * bt2 * mu2 - bt2 * bt2 * mu2 * dr.squared() + bt * bt2 * mu2 * secondDerivativeSumIdeal.getSum()
                    - Ans.getSumJEEMJEJE() + Ans.getSumUEE() - Ans.getSumJEMUExSquare() - Ans.getSumJEMUEySquare()
                    - Ans.getAEEJ0();
        }
        if (doIdeal) {
            x[1] = -nM * bt2 * mu2 - bt2 * bt2 * mu2 * dr.squared() + bt * bt2 * mu2 * secondDerivativeSumIdeal.getSum();
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
            vSumMinusIdeal.zeroSum();
            potentialMaster.calculate(box, allAtoms, vSumMinusIdeal);

            double bmu2 = bmu * bmu;
            for (int i = 0; i < nM; i++) {
                MoleculeAgent agentAtomI = leafAgentManager.getAgent(leafList.getAtom(i));
                IAtomOriented atom = (IAtomOriented) leafList.getAtom(i);

                dr.E(atom.getOrientation().getDirection());

                double ti = Math.atan2(dr.getX(1), dr.getX(0));
                double sinti = Math.sin(ti);
                double costi = Math.cos(ti);
                double cos2ti = Math.cos(2 * ti);

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

                agentAtomI.vEx().PE(vExiIdeal);
                agentAtomI.vEy().PE(vEyiIdeal);
                agentAtomI.vEEx().PE(vEExiIdeal);
                agentAtomI.vEEy().PE(vEEyiIdeal);
                agentAtomI.dvEx().PE(dvExidtiIdeal);
                agentAtomI.dvEy().PE(dvEyidtiIdeal);
                agentAtomI.dvEEx().PE(dvEExidtiIdeal);
                agentAtomI.dvEEy().PE(dvEEyidtiIdeal);
                agentAtomI.d2vEx().PE(d2vExidtidtiIdeal);
                agentAtomI.d2vEy().PE(d2vEyidtidtiIdeal);

            }

            vSumPairMinusIdeal.zeroSum();
            potentialMaster.calculate(box, allAtoms, vSumPairMinusIdeal);

            AEE = 0;
            JEMUEx = 0;
            JEMUEy = 0;
            JEEMJEJESelf = 0;
            UEESelf = 0;
            for (int i = 0; i < nM; i++) {
                MoleculeAgent agentAtomI = leafAgentManager.getAgent(leafList.getAtom(i));

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
                IAtomOriented atom = (IAtomOriented) leafList.getAtom(i);
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
        return new MoleculeAgent();
    }

    public void releaseAgent(MoleculeAgent agent, IAtom a, Box box) {

    }

}
