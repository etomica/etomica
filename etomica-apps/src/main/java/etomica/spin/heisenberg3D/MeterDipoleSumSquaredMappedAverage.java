package etomica.spin.heisenberg3D;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationEnergySum;
import etomica.potential.PotentialCalculationTorqueSum;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Null;

import java.io.Serializable;

/**
 * meter for AEE use mapping average
 *
 * @author Weisong
 */

public class MeterDipoleSumSquaredMappedAverage implements IEtomicaDataSource, AgentSource<MeterDipoleSumSquaredMappedAverage.MoleculeAgent> {

    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final Space space;
    protected final PotentialMaster potentialMaster;
    private final IteratorDirective allAtoms;
    //	private IBoundary boundary;
    protected PotentialCalculationEnergySum energySum;
    protected PotentialCalculationTorqueSum torqueSum;
    protected PotentialCalculationPhiSum secondDerivativeSum;
    protected Vector dr;
    protected Vector work;
    protected AtomLeafAgentManager atomAgentManager;
    protected AtomLeafAgentManager leafAgentManager;
    private Box box;
    private Vector vectorSum;
    //	private Vector r;
//	private Vector [] a;
    private double dipoleMagnitude;
    //	private double truncation;
    private double temperature;

    public MeterDipoleSumSquaredMappedAverage(final Space space, Box box, Simulation sim, double dipoleMagnitude, double temperature, PotentialMaster potentialMaster) {
        data = new DataDoubleArray(2);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{2});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.box = box;
        this.dipoleMagnitude = dipoleMagnitude;
        this.temperature = temperature;
        this.space = space;
        this.potentialMaster = potentialMaster;
        vectorSum = space.makeVector();
        vectorSum.setX(2, 1);

        leafAgentManager = new AtomLeafAgentManager<>(this, box, MoleculeAgent.class);
        torqueSum = new PotentialCalculationTorqueSum();
        energySum = new PotentialCalculationEnergySum();
        secondDerivativeSum = new PotentialCalculationPhiSum(space);
        torqueSum.setAgentManager(leafAgentManager);
        allAtoms = new IteratorDirective();
        dr = space.makeVector();
        work = space.makeVector();

        AtomLeafAgentManager.AgentSource<IntegratorVelocityVerlet.MyAgent> atomAgentSource = new AtomLeafAgentManager.AgentSource<IntegratorVelocityVerlet.MyAgent>() {
            public IntegratorVelocityVerlet.MyAgent makeAgent(IAtom a, Box agentBox) {
                return new IntegratorVelocityVerlet.MyAgent(space);
            }

            public void releaseAgent(IntegratorVelocityVerlet.MyAgent agent, IAtom atom, Box agentBox) {/**do nothing**/}
        };


    }

    public IData getData() {
//		IBoundary boundary = box.getBoundary();// TODO
        double[] x = data.getData();
        double bt = 1 / (temperature);//beta

        double mu = dipoleMagnitude;//miu
        double mu2 = mu * mu;
        double bt2 = bt * bt;
        double bt3 = bt * bt * bt;
        if (box == null) throw new IllegalStateException("no box");
        IAtomList leafList = box.getLeafList();

        int nM = leafList.getAtomCount();

        //TODO
//		IAtomOriented atom0 = (IAtomOriented) leafList.getAtom(0);
//		Vector  v0 =  (Vector) atom0.getOrientation().getDirection();
//		IAtomOriented atom1 = (IAtomOriented)  leafList.getAtom(0);
//		Vector  v1 =  (Vector) atom1.getOrientation().getDirection();
//		v1.E(0);
//		v1.setX(0, 1);
//		System.out.println("v0 = " + v0);
//		System.out.println("v1 = " + v1);

        torqueSum.reset();
        potentialMaster.calculate(box, allAtoms, torqueSum);
        secondDerivativeSum.zeroSum();
        potentialMaster.calculate(box, allAtoms, secondDerivativeSum);


        double A = 0;
        vectorSum.E(0);
        for (int i = 0; i < nM; i++) {
            IAtomOriented atom = (IAtomOriented) leafList.getAtom(0);
            dr.E(atom.getOrientation().getDirection());
            dr.normalize();
//			 A += -2.0/3.0*bt2*mu2*(dr.squared()-1);
            MoleculeAgent torqueAgent = (MoleculeAgent) leafAgentManager.getAgent(leafList.getAtom(i));
            dr.XE(torqueAgent.torque);
            vectorSum.PE(dr);
        }//i loop

        //x[0] = -0.25*bt2*bt2*mu2*vectorSum.squared(); //This part is zero
//        x[0] = -nM * bt2 * mu2 + 0.25 * bt3 * mu2 * secondDerivativeSum.getSum();
		x[0] = -nM*bt2*mu2 - 0.25*bt2*bt2*mu2*vectorSum.squared()+ 0.25*bt3*mu2*secondDerivativeSum.getSum();
        return data;
    }

//	public void setDipoleSource(DipoleSource newDipoleSource) {
//		dipoleSource = newDipoleSource;
//		secondDerivativeSum.setDipoleSource(newDipoleSource);
//	}

    public DataTag getTag() {
        return tag;
    }

    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    public Box getBox() {
        return box;
    }


    public MoleculeAgent makeAgent(IAtom a, Box box) {
        return new MoleculeAgent(space);
    }


    public void releaseAgent(MoleculeAgent agent, IAtom a, Box box) {

    }

    public Class getMoleculeAgentClass() {
        // TODO Auto-generated method stub
        return MoleculeAgent.class;
    }

    public AtomLeafAgentManager<IntegratorVelocityVerlet.MyAgent> getAtomAgentManager() {
        return atomAgentManager;
    }

    public static class MoleculeAgent implements Integrator.Torquable, Integrator.Forcible, Serializable {  //need public so to use with instanceof
        private static final long serialVersionUID = 1L;
        public final Vector torque;
        public final Vector force;

        public MoleculeAgent(Space space) {
            torque = space.makeVector();
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
