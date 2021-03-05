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
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationPhiSumHeisenberg;
import etomica.potential.PotentialCalculationTorqueSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterMappedAveraging implements IDataSource, AgentSource<MoleculeAgent> {
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final Space space;
    protected final PotentialMaster potentialMaster;
    protected final IteratorDirective allAtoms;
    //    protected PotentialCalculationFSum FSum;
    protected PotentialCalculationTorqueSum torqueSum;
    protected etomica.spin.heisenberg.PotentialCalculationPhiSum secondDerivativeSum;
    protected PotentialCalculationPhiSumHeisenberg secondDerivativeSumIdeal;
    protected double temperature;
    protected double J;
    protected double mu;
    protected double bt;
    protected int nMax;
    protected Vector dr, tmp;
    protected Vector work;
    protected AtomLeafAgentManager<MoleculeAgent> leafAgentManager;
    private Box box;
    protected PotentialCalculationHeisenberg Ans;

    public MeterMappedAveraging(Box box, double temperature, double interactionS, double dipoleMagnitude, PotentialMaster potentialMaster, int nMax) {
//        int a = 2*box.getLeafList().getAtomCount()+2;
        int nValues = 5;
        data = new DataDoubleArray(nValues);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{nValues});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.box = box;
        this.space = box.getSpace();
        this.temperature = temperature;
        this.potentialMaster = potentialMaster;
        J = interactionS;
        bt = 1 / temperature;
        mu = dipoleMagnitude;
        this.nMax = nMax;
        dr = space.makeVector();
        work = space.makeVector();
        tmp = space.makeVector();
        leafAgentManager = new AtomLeafAgentManager<MoleculeAgent>(this, box);
        torqueSum = new PotentialCalculationTorqueSum();
        torqueSum.setAgentManager(leafAgentManager);
//        FSum = new PotentialCalculationFSum(space, dipoleMagnitude, interactionS, bt);
        secondDerivativeSum = new PotentialCalculationPhiSum();
        secondDerivativeSum.setAgentManager(leafAgentManager);
        secondDerivativeSumIdeal = new PotentialCalculationPhiSumHeisenberg(space);

        this.nMax = nMax;
        Ans = new PotentialCalculationHeisenberg(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);
        allAtoms = new IteratorDirective();

    }

    public IData getData() {
        double[] x = data.getData();
        if (box == null) throw new IllegalStateException("no box");
        IAtomList leafList = box.getLeafList();
        torqueSum.reset();

//        MeterMappedAveraging.MoleculeAgent torqueAgent =  leafAgentManager.getAgent(leafList.getAtom(0));
//        double f1 = torqueAgent.torque.getX(1);
//        System.out.println("f1= "+f1);

        potentialMaster.calculate(box, allAtoms, torqueSum);

//        f1 = torqueAgent.torque.getX(1);
//        System.out.println("f1= "+f1);
//        System.exit(2);
//        FSum.zeroSum();
//        potentialMaster.calculate(box, allAtoms, FSum);
        secondDerivativeSum.reset();
        potentialMaster.calculate(box, allAtoms, secondDerivativeSum);

        Ans.zeroSum();
        potentialMaster.calculate(box, allAtoms, Ans);

        secondDerivativeSumIdeal.zeroSum();
        potentialMaster.calculate(box, allAtoms, secondDerivativeSumIdeal);

        double bt2 = bt * bt;
        double mu2 = mu * mu;
        int nM = leafList.size();
        double torqueScalar = 0;
        dr.E(0);
        tmp.E(0);
        for (int i = 0; i < nM; i++) {
            MoleculeAgent agentAtomI = leafAgentManager.getAgent(leafList.get(i));
            torqueScalar = agentAtomI.torque.getX(0);
            IAtomOriented atom = (IAtomOriented) leafList.get(i);
            dr.PEa1Tv1(torqueScalar, atom.getOrientation().getDirection());
            tmp.PE(atom.getOrientation().getDirection());
//            System.out.println("dr + " + dr);
        }//i loop
        x[0] = -nM * bt2 * mu2 - bt2 * bt2 * mu2 * dr.squared() + bt * bt2 * mu2 * secondDerivativeSumIdeal.getSum()
                - Ans.getSumJEEMJEJE() + Ans.getSumUEE() - Ans.getSumJEMUExSquare() - Ans.getSumJEMUEySquare()
                - Ans.getAEEJ0();
        x[1] = -nM * bt2 * mu2 - bt2 * bt2 * mu2 * dr.squared() + bt * bt2 * mu2 * secondDerivativeSumIdeal.getSum();
        x[2] = mu * tmp.getX(0);
        x[3] = mu * tmp.getX(1);
        x[4] = -Ans.getSumJEEMJEJE() + Ans.getSumUEE() - Ans.getSumJEMUExSquare() - Ans.getSumJEMUEySquare();
//        x[1] = Ans.getSumJEMUEx();
//        x[2] = Ans.getSumJEMUEy();
//        x[3] = Ans.getSumJEMUExIdeal();
//        x[4] = Ans.getSumJEMUEyIdeal();
//        x[5] = Ans.getSumJEEMJEJE();
//        x[6] = Ans.getSumUEE();
//        x[7] = -nM * bt2 * mu2 - bt2 * bt2 * mu2 * dr.squared() + bt * bt2 * mu2 * secondDerivativeSumIdeal.getSum();
//        x[8] = x[1] * x[1];
//        x[9] = x[2] * x[2];
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
