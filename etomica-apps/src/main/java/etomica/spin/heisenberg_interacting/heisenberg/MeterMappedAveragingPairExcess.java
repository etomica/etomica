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

/**
 * Computes pair-mapped average for a single spin pair
 */
public class MeterMappedAveragingPairExcess implements IDataSource, AgentSource<MoleculeAgent> {
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final Space space;
    protected final PotentialMaster potentialMaster;
    protected final IteratorDirective allAtoms;
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
    private final AtomPair pair;
    private final P2Spin p2Spin;
    protected final PotentialCalculationHeisenberg ans;

    public MeterMappedAveragingPairExcess(AtomPair pair, final Space space, Box box, Simulation sim, double temperature, double interactionS, double dipoleMagnitude, PotentialMaster potentialMaster, P2Spin p2Spin) {
//        int a = 2*box.getLeafList().getAtomCount()+2;
        this.pair = pair;
        this.p2Spin = p2Spin;
        int nValues = 8;
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

        dr = space.makeVector();
        work = space.makeVector();
        leafAgentManager = new AtomLeafAgentManager<MoleculeAgent>(this, box, MoleculeAgent.class);
        torqueSum = new PotentialCalculationTorqueSum();
        torqueSum.setAgentManager(leafAgentManager);
        secondDerivativeSum = new PotentialCalculationPhiSum();
        secondDerivativeSum.setAgentManager(leafAgentManager);
        secondDerivativeSumIdeal = new PotentialCalculationPhiSumHeisenberg(space);

        int nMax = 5;
        ans = new PotentialCalculationHeisenberg(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);
        allAtoms = new IteratorDirective();

    }

    public IData getData() {
        double[] x = data.getData();
        if (box == null) throw new IllegalStateException("no box");

        torqueSum.reset();
        potentialMaster.calculate(box, allAtoms, torqueSum);

        secondDerivativeSum.reset();
        secondDerivativeSum.doCalculation(pair,p2Spin);

        ans.zeroSum();
        ans.doCalculation(pair,p2Spin);

        secondDerivativeSumIdeal.zeroSum();
        secondDerivativeSumIdeal.doCalculation(pair,p2Spin);

        double bt2 = bt * bt;
        double mu2 = mu * mu;
        int nM = pair.getAtomCount();
        double torqueScalar = 0;
        dr.E(0);
        for (int i = 0; i < nM; i++) {
            MoleculeAgent agentAtomI = leafAgentManager.getAgent(pair.getAtom(i));
            torqueScalar = agentAtomI.torque.getX(0);
            IAtomOriented atom = (IAtomOriented) pair.getAtom(i);
            dr.PEa1Tv1(torqueScalar, atom.getOrientation().getDirection());
//            System.out.println("dr + " + dr);
        }//i loop
        x[0] = - ans.getSumJEEMJEJE() + ans.getSumUEE()
                - ans.getSumJEMUEx() * ans.getSumJEMUEx() - ans.getSumJEMUEy() * ans.getSumJEMUEy();
              //  - ans.getAEEJ0()  + ans.getSumJEMUExIdeal() * ans.getSumJEMUExIdeal() + ans.getSumJEMUEyIdeal() * ans.getSumJEMUEyIdeal();
        x[1] = ans.getSumJEMUEx();
        x[2] = ans.getSumJEMUEy();
        x[3] = ans.getSumJEMUExIdeal();//not used
        x[4] = ans.getSumJEMUEyIdeal();//not used
        x[5] = ans.getSumJEEMJEJE();//not used
//        x[6] = ans.getSumUEE();//not used//TODO
        x[6] =  ans.getAEEJ0();
//        x[0] -= -nM * bt2 * mu2 - bt2 * bt2 * mu2 * dr.squared() + bt * bt2 * mu2 * secondDerivativeSumIdeal.getSum();
        x[7] = -nM * bt2 * mu2 - bt2 * bt2 * mu2 * dr.squared() + bt * bt2 * mu2 * secondDerivativeSumIdeal.getSum();
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
