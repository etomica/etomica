package etomica.spin.heisenberg_interacting.heisenberg;

import etomica.atom.*;
import etomica.atom.AtomLeafAgentManager.AgentSource;
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

public class MeterMappedAveragingVSum3Pair implements IDataSource, AgentSource<MoleculeAgent> {
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
    protected PotentialCalculationMoleculeAgentSumPair vSumPair;
    protected double temperature;
    protected double J;
    protected double mu;
    protected double bt;
    protected Vector dr;
    protected Vector work;
    protected AtomLeafAgentManager<MoleculeAgent> leafAgentManager;
    private Box box;
    protected PotentialCalculationHeisenberg Ans;

    public MeterMappedAveragingVSum3Pair(final Space space, Box box, Simulation sim, double temperature, double interactionS, double dipoleMagnitude, IPotentialAtomic p2) {
//        int a = 2*box.getLeafList().getAtomCount()+2;
        int nValues = 24;
        data = new DataDoubleArray(nValues);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{nValues});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.box = box;
        this.p2 = p2;

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

        int nMax = 3;
        Ans = new PotentialCalculationHeisenberg(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);
        allAtoms = new IteratorDirective();
        vSum = new PotentialCalculationMoleculeAgentSum(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);
        vSumPair = new PotentialCalculationMoleculeAgentSumPair(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);

    }

    public IData getData() {
        double[] x = data.getData();
        if (box == null) throw new IllegalStateException("no box");
        IAtomList leafList = box.getLeafList();
        torqueSum.reset();
        secondDerivativeSum.reset();
        secondDerivativeSumIdeal.zeroSum();
//        MeterMappedAveraging.MoleculeAgent torqueAgent =  leafAgentManager.getAgent(leafList.getAtom(0));
//        double f1 = torqueAgent.torque.getX(1);
//        System.out.println("f1= "+f1);

        boolean twoPairOnly = true;

        AtomPair pair = new AtomPair();
        pair.atom0 = leafList.getAtom(0);
        pair.atom1 = leafList.getAtom(1);//01
        torqueSum.doCalculation(pair, p2);
        secondDerivativeSum.doCalculation(pair, p2);
        secondDerivativeSumIdeal.doCalculation(pair, p2);

        pair.atom1 = leafList.getAtom(2);//02
        torqueSum.doCalculation(pair, p2);
        secondDerivativeSum.doCalculation(pair, p2);
        secondDerivativeSumIdeal.doCalculation(pair, p2);

        if (!twoPairOnly) {
            pair.atom0 = leafList.getAtom(1);//12
            torqueSum.doCalculation(pair, p2);
            secondDerivativeSum.doCalculation(pair, p2);
            secondDerivativeSumIdeal.doCalculation(pair, p2);
        }

        Ans.zeroSum();
        pair.atom0 = leafList.getAtom(0);
        pair.atom1 = leafList.getAtom(1);//01
        Ans.doCalculation(pair, p2);

        pair.atom1 = leafList.getAtom(2);//02
        Ans.doCalculation(pair, p2);

        if (!twoPairOnly) {
            pair.atom0 = leafList.getAtom(1);//12
            Ans.doCalculation(pair, p2);
        }


        vSum.zeroSum();
        vSumPair.zeroSum();
        pair.atom0 = leafList.getAtom(0);
        pair.atom1 = leafList.getAtom(1);//01
        vSum.doCalculation(pair, p2);
        vSumPair.doCalculation(pair,p2);


        pair.atom1 = leafList.getAtom(2);//02
        vSum.doCalculation(pair, p2);
        vSumPair.doCalculation(pair,p2);

        if (!twoPairOnly) {
            pair.atom0 = leafList.getAtom(1);//12
            vSum.doCalculation(pair, p2);
            vSumPair.doCalculation(pair,p2);
        }


        double bt2 = bt * bt;
        double mu2 = mu * mu;
        double bmu = bt * mu;
        int nM = leafList.getAtomCount();


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
//            UEESelf -= vEEi * fi;
            //-fi*dvEi/dti*vEi
            double dvExi = agentAtomI.dvEx().getX(0);
            double dvEyi = agentAtomI.dvEy().getX(0);
//            UEESelf -= fi * (dvExi * vExi + dvEyi * vEyi);
            //vEi*phiii*vEi
            double phiii = bt * agentAtomI.phi.component(0, 0);
//            UEESelf += vExi * phiii * vExi + vEyi * phiii * vEyi;
            //-2*vEi*fEi
            //fExi = -bmu Sin[thetai]
            //fEyi = bmu Cos[thetai]
            IAtomOriented atom = (IAtomOriented) leafList.getAtom(i);
            double fExi = -bmu * atom.getOrientation().getDirection().getX(1);
            double fEyi = bmu * atom.getOrientation().getDirection().getX(0);
//            UEESelf -= 2 * (vExi * fExi + vEyi * fEyi);
            //-var[JEUME]
//            double tmp = dvExi + bmu * atom.getOrientation().getDirection().getX(0) + vExi * fi;
//            JEMUEx += tmp*tmp;
//            tmp =dvEyi + bmu * atom.getOrientation().getDirection().getX(1) + vEyi * fi;
//            JEMUEy += tmp*tmp;

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
        }//i loop
        x[0] = -nM * bt2 * mu2 - bt2 * bt2 * mu2 * dr.squared() + bt * bt2 * mu2 * secondDerivativeSumIdeal.getSum()
                - Ans.getSumJEEMJEJE() + Ans.getSumUEE() - Ans.getSumJEMUExSquare() - Ans.getSumJEMUEySquare()
                - Ans.getAEEJ0();

        x[1] = Ans.getSumJEMUEx();
        x[2] = Ans.getSumJEMUEy();
        x[3] = Ans.getSumJEMUExIdeal();
        x[4] = Ans.getSumJEMUEyIdeal();
        x[5] = Ans.getSumJEEMJEJE();
        x[6] = Ans.getSumUEE();
        x[7] = -nM * bt2 * mu2 - bt2 * bt2 * mu2 * dr.squared() + bt * bt2 * mu2 * secondDerivativeSumIdeal.getSum();
        x[8] = Ans.getSumJEMUEx() * Ans.getSumJEMUEx();
        x[9] = Ans.getSumJEMUEy() * Ans.getSumJEMUEy();
        x[10] = -x[5] + x[6] - Ans.getSumJEMUExSquare() - Ans.getSumJEMUEySquare();
        x[11] = Ans.getSumJEMUExIdeal() * Ans.getSumJEMUExIdeal();
        x[12] = Ans.getSumJEMUEyIdeal() * Ans.getSumJEMUEyIdeal();
        x[13] = Ans.getAEEJ0();
        x[14] = -nM * bt2 * mu2;
        x[15] = -bt2 * bt2 * mu2 * dr.squared();
        x[16] = bt * bt2 * mu2 * secondDerivativeSumIdeal.getSum();


        x[17] = AEE;
        x[18] = JEMUEx;
        x[19] = JEMUEy;
        x[20] = JEEMJEJESelf + vSumPair.getSumJEEMJEJE();
//        x[21] = UEESelf + vSum.getSumUEE();
        x[21] = vSumPair.getSumUEE() ;
        x[22] = JEMUEx * JEMUEx;
        x[23] = JEMUEy * JEMUEy;
        return data;
    }

    //    int count = 0;
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
