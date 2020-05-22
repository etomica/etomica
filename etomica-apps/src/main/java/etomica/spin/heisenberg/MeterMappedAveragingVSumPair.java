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
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterMappedAveragingVSumPair implements IDataSource, AgentSource<MoleculeAgent> {
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final Space space;
    protected final IPotentialAtomic p2;
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
    protected PotentialCalculationPair pair;
    protected double temperature;
    protected double J;
    protected double mu;
    protected double bt;
    protected int nMax;
    protected Vector dr;
    protected Vector work;
    protected AtomLeafAgentManager<MoleculeAgent> leafAgentManager;
    private Box box;
    protected PotentialCalculationHeisenberg Ans;

    public MeterMappedAveragingVSumPair(final Space space, Box box, Simulation sim, double temperature, double interactionS, double dipoleMagnitude, IPotentialAtomic p2, int nMax) {
//        int a = 2*box.getLeafList().getAtomCount()+2;
        int nValues = 24;
        data = new DataDoubleArray(nValues);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{nValues});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.box = box;
        this.p2 = p2;
        this.nMax = nMax;
        this.space = space;
        this.temperature = temperature;
        J = interactionS;
        bt = 1 / temperature;
        mu = dipoleMagnitude;

        dr = space.makeVector();
        work = space.makeVector();
        leafAgentManager = new AtomLeafAgentManager<MoleculeAgent>(this, box, MoleculeAgent.class);
        torqueSum = new PotentialCalculationTorqueSum();
        torqueSum.setAgentManager(leafAgentManager);
//        FSum = new PotentialCalculationFSum(space, dipoleMagnitude, interactionS, bt);
        energySum = new PotentialCalculationEnergySum();
        secondDerivativeSum = new PotentialCalculationPhiSum();
        secondDerivativeSum.setAgentManager(leafAgentManager);

        secondDerivativeSumIdeal = new PotentialCalculationPhiSumHeisenberg(space);


        Ans = new PotentialCalculationHeisenberg(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);
        vSum = new PotentialCalculationMoleculeAgentSum(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);
        vSumPair = new PotentialCalculationMoleculeAgentSumPair(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);

        vSumMinusIdeal = new PotentialCalculationMoleculeAgentSumMinusIdeal(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);
        vSumPairMinusIdeal = new PotentialCalculationMoleculeAgentSumMinusIdealPair(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);


        pair = new PotentialCalculationPair(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);

        allAtoms = new IteratorDirective();

    }

    public IData getData() {
        double[] x = data.getData();
        if (box == null) throw new IllegalStateException("no box");
        IAtomList leafList = box.getLeafList();
        double bt2 = bt * bt;
        double mu2 = mu * mu;
        double bmu = bt * mu;
        int nM = leafList.getAtomCount();

        torqueSum.reset();
        torqueSum.doCalculation(box.getLeafList(), p2);
//        f1 = torqueAgent.torque.getX(1);
//        System.out.println("f1= "+f1);
//        System.exit(2);
//        FSum.zeroSum();
//        potentialMaster.calculate(box, allAtoms, FSum);
        secondDerivativeSum.reset();
        secondDerivativeSum.doCalculation(box.getLeafList(), p2);


        Ans.zeroSum();
        Ans.doCalculation(box.getLeafList(), p2);

        vSum.zeroSum();
        vSum.doCalculation(box.getLeafList(), p2);
        vSumPair.zeroSum();
        vSumPair.doCalculation(box.getLeafList(), p2);

        pair.zeroSum();
        pair.doCalculation(box.getLeafList(), p2);

        secondDerivativeSumIdeal.zeroSum();
        secondDerivativeSumIdeal.doCalculation(box.getLeafList(), p2);


        double AEE = 0, JEMUEx = 0, JEMUEy = 0;
        double JEEMJEJESelf = 0, UEESelf = 0;
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


        AEE = -vSumPair.getSumJEEMJEJE() - JEEMJEJESelf + UEESelf + vSumPair.getSumUEE() - JEMUEx * JEMUEx - JEMUEy * JEMUEy;


        double torqueScalar = 0;
        dr.E(0);
        for (int i = 0; i < nM; i++) {
            MoleculeAgent agentAtomI = leafAgentManager.getAgent(leafList.getAtom(i));
            torqueScalar = agentAtomI.torque.getX(0);
            IAtomOriented atom = (IAtomOriented) leafList.getAtom(i);
            dr.PEa1Tv1(torqueScalar, atom.getOrientation().getDirection());
//            System.out.println("dr + " + dr);
        }//i loop
        x[0] = -nM * bt2 * mu2 - bt2 * bt2 * mu2 * dr.squared() + bt * bt2 * mu2 * secondDerivativeSumIdeal.getSum()
                - Ans.getSumJEEMJEJE() + Ans.getSumUEE() - Ans.getSumJEMUEx() * Ans.getSumJEMUEx() - Ans.getSumJEMUEy() * Ans.getSumJEMUEy()
                - Ans.getAEEJ0();

        x[1] = Ans.getSumJEMUEx();
        x[2] = Ans.getSumJEMUEy();
        x[3] = Ans.getSumJEMUExIdeal();
        x[4] = Ans.getSumJEMUEyIdeal();

        x[5] = Ans.getSumJEEMJEJE();
        x[6] = Ans.getSumUEE();
        x[7] = -nM * bt2 * mu2 - bt2 * bt2 * mu2 * dr.squared() + bt * bt2 * mu2 * secondDerivativeSumIdeal.getSum();
        x[8] = x[1] * x[1];
        x[9] = x[2] * x[2];


        x[10] = AEE;
        x[11] = JEMUEx;
        x[12] = JEMUEy;
        x[13] = JEEMJEJESelf + vSumPair.getSumJEEMJEJE();
        x[14] = UEESelf + vSumPair.getSumUEE();
        x[15] = JEMUEx * JEMUEx;
        x[16] = JEMUEy * JEMUEy;


        vSumMinusIdeal.zeroSum();
        vSumMinusIdeal.doCalculation(box.getLeafList(), p2);
        double bmu2 = bmu * bmu;
        for (int i = 0; i < nM; i++) {
            MoleculeAgent agentAtomI = leafAgentManager.getAgent(leafList.getAtom(i));
            IAtomOriented atom = (IAtomOriented) leafList.getAtom(i);

            dr.E(atom.getOrientation().getDirection());

            double ti = Math.atan2(dr.getX(1), dr.getX(0));
            double sinti = Math.sin(ti);
            double costi = Math.cos(ti);
            double cos2ti = Math.cos(2 * ti);
            double sin2ti = Math.sin(2 * ti);

            double vExiIdeal = -bmu * sinti;
            double vEyiIdeal = bmu * costi;
            double vEExiIdeal = 0.5 * bmu2 * costi * sinti;
            double vEEyiIdeal = -0.5 * bmu2 * costi * sinti;
            double dvExidtiIdeal = -bmu * costi;
            double dvEyidtiIdeal = -bmu * sinti;
            double dvEExidtiIdeal = 0.5 * bmu2 * cos2ti;
            double dvEEyidtiIdeal = -dvEExidtiIdeal;
            double d2vExidtidtiIdeal = bmu * sinti;
            double d2vEyidtidtiIdeal = -bmu*costi;

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
        vSumPairMinusIdeal.doCalculation(box.getLeafList(), p2);

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


        x[17] = AEE;
        x[18] = JEMUEx;
        x[19] = JEMUEy;
        x[20] = JEEMJEJESelf + vSumPairMinusIdeal.getSumJEEMJEJE();
        x[21] = UEESelf + vSumPairMinusIdeal.getSumUEE();
        x[22] = JEMUEx * JEMUEx;
        x[23] = JEMUEy * JEMUEy;


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


}
