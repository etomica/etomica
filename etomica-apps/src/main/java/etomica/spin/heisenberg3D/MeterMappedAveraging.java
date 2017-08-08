package etomica.spin.heisenberg3D;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.Integrator;
import etomica.molecule.DipoleSource;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space1d.Vector1D;
import etomica.units.Null;

import java.io.Serializable;

public class MeterMappedAveraging implements IEtomicaDataSource, AgentSource<MeterMappedAveraging.MoleculeAgent> {


    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final Space space;
    protected final PotentialMaster potentialMaster;
    protected final IteratorDirective allAtoms;
    public double QValue = 0;
    //private IBoundary boundary;
    protected PotentialCalculationEnergySum energySum;
    protected PotentialCalculationFSum FSum;
    protected PotentialCalculationTorqueSum torqueSum;
    protected PotentialCalculationPhiSumHeisenberg secondDerivativeSum;
    //private double truncation;
    protected double temperature;
    protected double J;
    protected double mu;
    protected double bt;
    protected Vector dr;
    protected Vector work;
    protected AtomLeafAgentManager leafAgentManager;
    protected DipoleSource dipoleSource;
    private Box box;


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

        dr = new Vector1D();
        leafAgentManager = new AtomLeafAgentManager<MoleculeAgent>(this, box, MoleculeAgent.class);
        torqueSum = new PotentialCalculationTorqueSum();
        torqueSum.setAgentManager(leafAgentManager);
        FSum = new PotentialCalculationFSum(space, dipoleMagnitude, interactionS, temperature);
        energySum = new PotentialCalculationEnergySum();
        secondDerivativeSum = new PotentialCalculationPhiSumHeisenberg(space);


        allAtoms = new IteratorDirective();

    }

    @Override
    public IData getData() {
        double[] x = data.getData();
        double bt = 1 / (temperature);//beta
        if (box == null) throw new IllegalStateException("no box");

        IAtomList leafList = box.getLeafList();
        torqueSum.reset();
        potentialMaster.calculate(box, allAtoms, torqueSum);
        secondDerivativeSum.zeroSum();
        potentialMaster.calculate(box, allAtoms, secondDerivativeSum);


        double bt2 = bt * bt;
        double mu2 = mu * mu;
        int nM = leafList.getAtomCount();
        double A = 0;
        dr.E(0);

        //torque square sum is zero so don't need it here.
		 for (int i = 0;i < nM; i++){
			 MoleculeAgent torqueAgent = (MoleculeAgent) leafAgentManager.getAgent(leafList.getAtom(i));
			 dr.PE(torqueAgent.torque);
//            System.out.println(torqueAgent.torque);

			 //test for <f(1-x^2)> the result is zero!!!!!
//			 IAtomOriented atom = (IAtomOriented)leafList.getAtom(0);
//			 double ex = atom.getOrientation().getDirection().getX(0);
//			 double ey = atom.getOrientation().getDirection().getX(1);
//			 torqueSum.PEa1Tv1((ex+ey),torqueAgent.torque);
		 }//i loop

		x[0] = -nM*bt2*mu2 - 0.25*bt2*bt2*mu2*dr.squared()+ 0.25*J*bt*bt2*mu2*secondDerivativeSum.getSum();
//        x[0] = -nM * bt2 * mu2 + 0.25 * J * bt * bt2 * mu2 * secondDerivativeSum.getSum();

//		test for <f(1-x^2)>  the result is zero!!!!!
//		x[0] = dr.squared();
//		x[0] = secondDerivativeSum.getSum();
        return data;
    }


    public DataTag getTag() {
        return tag;
    }

    public IEtomicaDataInfo getDataInfo() {
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

    public static class MoleculeAgent implements Integrator.Torquable, Integrator.Forcible, Serializable {  //need public so to use with instanceof
        private static final long serialVersionUID = 1L;
        public final Vector torque;
        public final Vector force;

        public MoleculeAgent(Space space) {
            torque = new Vector1D();
            force = space.makeVector();
        }

        public Vector torque() {
            return torque;
        }

        public Vector force() {
            return force;
        }
    }

}
