package etomica.spin.heisenberg_interacting.heisenberg;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.Integrator;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space1d.Tensor1D;
import etomica.space1d.Vector1D;
import etomica.space2d.Vector2D;
import etomica.units.dimensions.Null;

public class MeterMappedAveragingPair implements IDataSource, AgentSource<MeterMappedAveraging.MoleculeAgent> {
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
    protected double temperature;
    protected double J;
    protected double mu;
    protected double bt;
    protected Vector dr;
    protected Vector work;
    protected AtomLeafAgentManager<MeterMappedAveraging.MoleculeAgent> leafAgentManager;
    private Box box;
    protected PotentialCalculationHeisenberg Ans;

    public MeterMappedAveragingPair(final Space space, Box box, Simulation sim, double temperature, double interactionS, double dipoleMagnitude, IPotentialAtomic p2) {
//        int a = 2*box.getLeafList().getAtomCount()+2;
        data = new DataDoubleArray(6);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{6});
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
        leafAgentManager = new AtomLeafAgentManager<MeterMappedAveraging.MoleculeAgent>(this, box, MeterMappedAveraging.MoleculeAgent.class);
        torqueSum = new PotentialCalculationTorqueSum();
        torqueSum.setAgentManager(leafAgentManager);
//        FSum = new PotentialCalculationFSum(space, dipoleMagnitude, interactionS, bt);
        energySum = new PotentialCalculationEnergySum();
        secondDerivativeSum = new PotentialCalculationPhiSum();
        secondDerivativeSum.setAgentManager(leafAgentManager);


        int nMax = 50;
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

        torqueSum.doCalculation(box.getLeafList(),p2);

//        f1 = torqueAgent.torque.getX(1);
//        System.out.println("f1= "+f1);
//        System.exit(2);
//        FSum.zeroSum();
//        potentialMaster.calculate(box, allAtoms, FSum);
        secondDerivativeSum.reset();
        secondDerivativeSum.doCalculation(box.getLeafList(),p2);

        Ans.zeroSum();
        Ans.doCalculation(box.getLeafList(),p2);

//        x[1] = Ans.getSumJEMUEx() * Ans.getSumJEMUEx();
//        x[2] = Ans.getSumJEMUEy() * Ans.getSumJEMUEy();
//        x[3] = Ans.getSumJEMUEx();
//        x[4] = Ans.getSumJEMUEy();
        x[0] = Ans.getSumJEEMJEJE();
        x[1] = Ans.getSumUEE();
        x[2] = Ans.getdipolex() + Ans.getdipoley();
        x[3] = Ans.getSumJEMUE() * Ans.getSumJEMUE();
        x[4] = Ans.getSumJEMUE();
        x[5] = Ans.getdipoleconv();

        return data;

    }


    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public MeterMappedAveraging.MoleculeAgent makeAgent(IAtom a, Box box) {
        return new MeterMappedAveraging.MoleculeAgent();
    }

    public void releaseAgent(MeterMappedAveraging.MoleculeAgent agent, IAtom a, Box box) {

    }



}
