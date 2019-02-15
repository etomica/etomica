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

public class MeterMappedAveraging3Pair implements IDataSource, AgentSource<MoleculeAgent> {
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
    protected double temperature;
    protected double J;
    protected double mu;
    protected double bt;
    protected Vector dr;
    protected Vector work;
    protected AtomLeafAgentManager<MoleculeAgent> leafAgentManager;
    private Box box;
    protected PotentialCalculationHeisenberg Ans;

    public MeterMappedAveraging3Pair(final Space space, Box box, Simulation sim, double temperature, double interactionS, double dipoleMagnitude, IPotentialAtomic p2) {
//        int a = 2*box.getLeafList().getAtomCount()+2;
        int nValues = 17;
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

        if(!twoPairOnly) {
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

        if(!twoPairOnly) {
            pair.atom0 = leafList.getAtom(1);//12
            Ans.doCalculation(pair, p2);
        }

        double bt2 = bt * bt;
        double mu2 = mu * mu;
        int nM = leafList.getAtomCount();
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
        x[15] = - bt2 * bt2 * mu2 * dr.squared();
        x[16] = bt * bt2 * mu2 * secondDerivativeSumIdeal.getSum();
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
