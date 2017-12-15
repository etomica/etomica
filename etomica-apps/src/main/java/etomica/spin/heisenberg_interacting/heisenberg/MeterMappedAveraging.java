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
import etomica.units.dimensions.Null;

import java.io.Serializable;

public class MeterMappedAveraging implements IDataSource, AgentSource<MeterMappedAveraging.MoleculeAgent> {


    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final Space space;
    protected final PotentialMaster potentialMaster;
    protected final IteratorDirective allAtoms;
    //private IBoundary boundary;
    protected PotentialCalculationEnergySum energySum;
    //    protected PotentialCalculationFSum FSum;
//    protected PotentialCalculationTorqueSum torqueSum;
//    protected PotentialCalculationPhiSumHeisenberg secondDerivativeSum;
    //private double truncation;
    protected double temperature;
    protected double J;
    protected double mu;
    protected double bt;
    protected Vector dr;
    protected Vector work;
    protected AtomLeafAgentManager<MoleculeAgent> leafAgentManager;
    private Box box;

    //TODO debug only
    protected PotentialCalculationHeisenberg Ans;

    public MeterMappedAveraging(final Space space, Box box, Simulation sim, double temperature, double interactionS, double dipoleMagnitude, PotentialMaster potentialMaster) {
        data = new DataDoubleArray(2);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{2});
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
//        torqueSum = new PotentialCalculationTorqueSum();
//        torqueSum.setAgentManager(leafAgentManager);
//        FSum = new PotentialCalculationFSum(space, dipoleMagnitude, interactionS, bt);
        energySum = new PotentialCalculationEnergySum();
//        secondDerivativeSum = new PotentialCalculationPhiSumHeisenberg(space, J, bt);

        int nMax = 3;
        Ans = new PotentialCalculationHeisenberg(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);
        allAtoms = new IteratorDirective();

    }

    public IData getData() {
        double[] x = data.getData();
        if (box == null) throw new IllegalStateException("no box");
        IAtomList leafList = box.getLeafList();
//        torqueSum.reset();
//        potentialMaster.calculate(box, allAtoms, torqueSum);
//        FSum.zeroSum();
//        potentialMaster.calculate(box, allAtoms, FSum);
//        secondDerivativeSum.zeroSum();
//        potentialMaster.calculate(box, allAtoms, secondDerivativeSum);

        Ans.zeroSumJEEMJEJE();
        Ans.zeroSumJEMUE();
        Ans.zeroSumJEMUESquare();
        Ans.zeroSumUEE();
        potentialMaster.calculate(box, allAtoms, Ans);
        int nM = leafList.getAtomCount();

//        x[0] = (-Ans.getSumJEEMJEJE() + Ans.getSumUEE() - Ans.getSumJEMUESquare())+3*nM*bt*bt;
        x[0] = (-Ans.getSumJEEMJEJE() + Ans.getSumUEE() - Ans.getSumJEMUESquare());
        x[1] = Ans.getSumJEMUE();
        return data;
    }


    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public MoleculeAgent makeAgent(IAtom a, Box box) {
        return new MoleculeAgent(space);
    }

    public void releaseAgent(MoleculeAgent agent, IAtom a, Box box) {

    }

    public Class getAgentClass() {
        return MoleculeAgent.class;
    }

    public static class MoleculeAgent {  //need public so to use with instanceof
        public final Vector torque;
        public final Tensor phi;

        public MoleculeAgent(Space space) {
            torque = new Vector1D();
            phi = new Tensor1D();
        }

        public Vector torque() {
            return torque;
        }

        public Tensor phi() {
            return phi;
        }
    }

}
