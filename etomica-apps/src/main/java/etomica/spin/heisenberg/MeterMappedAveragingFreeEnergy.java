package etomica.spin.heisenberg;

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
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterMappedAveragingFreeEnergy implements IDataSource, AgentSource<MoleculeAgent> {
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final Space space;
    protected final PotentialMaster potentialMaster;
    protected final IteratorDirective allAtoms;
    protected PotentialCalculationEnergySum energySum;
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
    protected PotentialCalculationFreeEnergy FE;

    public MeterMappedAveragingFreeEnergy(final Space space, Box box, Simulation sim, double temperature, double interactionS, double dipoleMagnitude, PotentialMaster potentialMaster) {
        int nValues = 5;
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
        tmp = space.makeVector();
        leafAgentManager = new AtomLeafAgentManager<MoleculeAgent>(this, box, MoleculeAgent.class);
        torqueSum = new PotentialCalculationTorqueSum();
        torqueSum.setAgentManager(leafAgentManager);
        energySum = new PotentialCalculationEnergySum();
        secondDerivativeSum = new PotentialCalculationPhiSum();
        secondDerivativeSum.setAgentManager(leafAgentManager);
        secondDerivativeSumIdeal = new PotentialCalculationPhiSumHeisenberg(space);

        FE = new PotentialCalculationFreeEnergy(space, dipoleMagnitude, interactionS, bt, leafAgentManager);
        allAtoms = new IteratorDirective();

    }

    public IData getData() {
        double[] x = data.getData();
        if (box == null) throw new IllegalStateException("no box");
        IAtomList leafList = box.getLeafList();
        torqueSum.reset();
        potentialMaster.calculate(box, allAtoms, torqueSum);


        FE.zeroSum();
        potentialMaster.calculate(box, allAtoms, FE);
        x[0] =  FE.getSumU_Map();

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
