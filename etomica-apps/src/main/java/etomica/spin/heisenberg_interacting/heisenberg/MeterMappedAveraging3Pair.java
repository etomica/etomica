package etomica.spin.heisenberg_interacting.heisenberg;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.AtomPair;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.potential.IPotentialAtomic;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationEnergySum;
import etomica.potential.PotentialCalculationTorqueSum;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterMappedAveraging3Pair implements IDataSource, AgentSource<MeterMappedAveraging.MoleculeAgent> {
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

    public MeterMappedAveraging3Pair(final Space space, Box box, Simulation sim, double temperature, double interactionS, double dipoleMagnitude, IPotentialAtomic p2) {
//        int a = 2*box.getLeafList().getAtomCount()+2;
        data = new DataDoubleArray(2);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{2});
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


        int nMax = 10;
        Ans = new PotentialCalculationHeisenberg(space, dipoleMagnitude, interactionS, bt, nMax, leafAgentManager);
        allAtoms = new IteratorDirective();

    }

    public IData getData() {
        double[] x = data.getData();
        if (box == null) throw new IllegalStateException("no box");
        IAtomList leafList = box.getLeafList();
        torqueSum.reset();
        secondDerivativeSum.reset();
//        MeterMappedAveraging.MoleculeAgent torqueAgent =  leafAgentManager.getAgent(leafList.getAtom(0));
//        double f1 = torqueAgent.torque.getX(1);
//        System.out.println("f1= "+f1);

        AtomPair pair = new AtomPair();
        pair.atom0 = box.getLeafList().getAtom(0);
        pair.atom1 = box.getLeafList().getAtom(1);//01

        torqueSum.doCalculation(pair,p2);
        secondDerivativeSum.doCalculation(pair,p2);
        pair.atom1 = box.getLeafList().getAtom(2);//02
        torqueSum.doCalculation(pair,p2);
        secondDerivativeSum.doCalculation(pair,p2);
        pair.atom0 = box.getLeafList().getAtom(1);//12
        torqueSum.doCalculation(pair,p2);//12
        secondDerivativeSum.doCalculation(pair,p2);

        Ans.zeroSum();
        pair.atom0 = box.getLeafList().getAtom(0);
        pair.atom1 = box.getLeafList().getAtom(1);//01
        Ans.doCalculation(pair,p2);
        pair.atom1 = box.getLeafList().getAtom(2);//02
        Ans.doCalculation(pair,p2);
        pair.atom0 = box.getLeafList().getAtom(1);//12
        Ans.doCalculation(pair,p2);


        x[0] = - Ans.getSumJEEMJEJE() + Ans.getSumUEE() - Ans.getSumJEMUE() * Ans.getSumJEMUE();
        x[1]  = Ans.getSumJEMUE();
//        System.out.println("mine "+x[0]+" "+x[1]);


//count += 1;
//if(count > 19){
//    System.exit(2);
//}
        return data;

    }

//    int count = 0;
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
