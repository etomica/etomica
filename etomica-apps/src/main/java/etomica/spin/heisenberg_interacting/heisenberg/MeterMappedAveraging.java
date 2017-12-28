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

public class MeterMappedAveraging implements IDataSource, AgentSource<MeterMappedAveraging.MoleculeAgent> {


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
    protected double temperature;
    protected double J;
    protected double mu;
    protected double bt;
    protected Vector dr;
    protected Vector work;
    protected AtomLeafAgentManager<MoleculeAgent> leafAgentManager;
    private Box box;
    protected PotentialCalculationHeisenberg Ans;

    public MeterMappedAveraging(final Space space, Box box, Simulation sim, double temperature, double interactionS, double dipoleMagnitude, PotentialMaster potentialMaster) {
        data = new DataDoubleArray(5);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{5});
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
//        FSum = new PotentialCalculationFSum(space, dipoleMagnitude, interactionS, bt);
        energySum = new PotentialCalculationEnergySum();
        secondDerivativeSum = new PotentialCalculationPhiSum();
        secondDerivativeSum.setAgentManager(leafAgentManager);


        int nMax = 10;
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

        Ans.zeroSumJEEMJEJE();
        Ans.zeroSumJEMUEx();
        Ans.zeroSumJEMUEy();
        Ans.zeroSumJEMUESquare();
        Ans.zeroSumUEE();

//        System.out.println(Ans.getSumJEMUEx());
//        System.out.println(Ans.getSumJEMUEy());
//        System.exit(2);

        potentialMaster.calculate(box, allAtoms, Ans);

//        int nM = leafList.getAtomCount();
//        x[0] = (-Ans.getSumJEEMJEJE() + Ans.getSumUEE() - Ans.getSumJEMUESquare())+3*nM*bt*bt;

        x[0] = Ans.getSumJEMUEx();
        x[1] = Ans.getSumJEMUEy();
        x[2] = x[0] * x[0];
        x[3] = x[1] * x[1];
        x[4] = -Ans.getSumJEEMJEJE() + Ans.getSumUEE();

        //TODO debug only
//        x[0] = Ans.getSumJEMUEx();
//        x[1] = Ans.getSumJEMUEy();
//        x[2] = Ans.getSumJEMUESquare();
//        x[3] = 0;
//        x[4] = -Ans.getSumJEEMJEJE() + Ans.getSumUEE();
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


    public static class MoleculeAgent implements Integrator.Torquable, Integrator.Forcible {  //need public so to use with instanceof
        public final Vector torque;
        public final Tensor phi;
        public final Vector force;

        public MoleculeAgent() {
            torque = new Vector1D();
            phi = new Tensor1D();
            force = new Vector2D();
        }

        public Vector torque() {
            return torque;
        }

        public Tensor phi() {
            return phi;
        }

        public Vector force() {
            return force;
        }
    }

}
