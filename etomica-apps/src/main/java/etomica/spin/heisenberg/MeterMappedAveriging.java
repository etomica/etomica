package etomica.spin.heisenberg;

import java.io.Serializable;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.ISimulation;
import etomica.api.IVectorMutable;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.DipoleSource;
import etomica.atom.IAtomOriented;
import etomica.atom.iterator.IteratorDirective;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.PotentialCalculationEnergySum;
import etomica.potential.PotentialCalculationFSum;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialCalculationPhiSumHeisenberg;
import etomica.space.ISpace;
import etomica.space1d.Vector1D;
import etomica.units.Null;
import etomica.util.numerical.BesselFunction;

public class MeterMappedAveriging implements IEtomicaDataSource ,AgentSource<MeterMappedAveriging.MoleculeAgent> {
	

	protected final DataDoubleArray data;
	protected final DataInfoDoubleArray dataInfo;
	protected final DataTag tag;
	//private IBoundary boundary;
	protected PotentialCalculationEnergySum energySum;
	protected PotentialCalculationFSum FSum;
	protected PotentialCalculationPhiSumHeisenberg secondDerivativeSum;
	protected final ISpace space;
	private IBox box;
	private IVectorMutable torqueSum;
	//private IVectorMutable r;
	//private IVectorMutable [] a;
	
	//private double truncation;
	protected double temperature;
	protected double J;
	protected double mu;
	protected double bt;

	public double QValue = 0;
	protected final IPotentialMaster potentialMaster;
	protected final IteratorDirective allAtoms;
	protected IVectorMutable dr;
	protected IVectorMutable work;
	protected AtomLeafAgentManager leafAgentManager;
	protected DipoleSource dipoleSource;
	protected AtomLeafAgentManager atomAgentManager;
	protected PotentialCalculationForceSum pcForce;

	
	//use torquesum here TODO
	public MeterMappedAveriging(final ISpace space, IBox box, ISimulation sim, double temperature,double interactionS,double dipoleMagnitude,IPotentialMaster potentialMaster) {
        data = new DataDoubleArray(2);
		dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{2}); 
		tag = new DataTag();
		dataInfo.addTag(tag);
		this.box = box;
		this.space = space;
		this.temperature = temperature;
		this.potentialMaster = potentialMaster;
		J = interactionS;
		bt = 1/temperature;
		mu = dipoleMagnitude;
//        QValue = bt*bt*mu*mu + org.apache.commons.math3.b ;//How use Bessel fucntion here????
		
		torqueSum = space.makeVector();
		torqueSum.setX(2, 1);
		FSum = new PotentialCalculationFSum(space,dipoleMagnitude,interactionS,temperature);
		energySum = new PotentialCalculationEnergySum();
		secondDerivativeSum = new  PotentialCalculationPhiSumHeisenberg(space,dipoleMagnitude,interactionS,temperature);
		leafAgentManager  = new AtomLeafAgentManager<MoleculeAgent>(this , box, MoleculeAgent.class); 
//		FSum.setAgentManager(leafAgentManager); TODO
		
		allAtoms = new IteratorDirective();
		dr = space.makeVector();
		work = space.makeVector();
		
		AtomLeafAgentManager.AgentSource<IntegratorVelocityVerlet.MyAgent> atomAgentSource = new AtomLeafAgentManager.AgentSource<IntegratorVelocityVerlet.MyAgent>() {
		    public IntegratorVelocityVerlet.MyAgent makeAgent(IAtom a,IBox box) {
		        return new IntegratorVelocityVerlet.MyAgent(space);
		    }
		    public void releaseAgent(IntegratorVelocityVerlet.MyAgent agent, IAtom atom,IBox box) {/**do nothing**/}
        };
	}

	@Override
	public IData getData() {
		double[] x = data.getData();
		double bt = 1/(temperature);//beta
		if (box == null) throw new IllegalStateException("no box");
		
		IAtomList leafList = box.getLeafList();
		 FSum.zeroSum();
		 potentialMaster.calculate(box, allAtoms, FSum);
		 secondDerivativeSum.zeroSum();
		 potentialMaster.calculate(box, allAtoms, secondDerivativeSum);
		 
		
		 double bt2 = bt*bt;
		 double bt3 = bt2*bt;
		 double mu2 = mu*mu;
		 double mu3 = mu2*mu;
		 
		 
		 double Q = bt*bt*mu*mu*(1+BesselFunction.I(1, J*bt)/BesselFunction.I(0, J*bt));
		 IAtomOriented atom1 = (IAtomOriented) leafList.getAtom(0);
		 IAtomOriented atom2 = (IAtomOriented) leafList.getAtom(1);
		 
		 
		 int nM = leafList.getAtomCount();	
		 double A = 0;
		 torqueSum.E(0);
		 for (int i = 0;i < nM; i++){

			 IAtomOriented atom = (IAtomOriented) leafList.getAtom(i);
			 MoleculeAgent torqueAgent = (MoleculeAgent) leafAgentManager.getAgent(leafList.getAtom(i));
			 torqueSum.PE(torqueAgent.torque);
		 }//i loop


		x[0] = -nM*bt2*mu2 - 0.25*bt2*bt2*mu2*torqueSum.squared()+ secondDerivativeSum.getSum();
//		 x[0] = secondDerivativeSum.getSum();//TODO
//		 x[1] = 1;
		
		return data;
	}
	
	public void setDipoleSource(DipoleSource newDipoleSource) {
		dipoleSource = newDipoleSource;
		secondDerivativeSum.setDipoleSource(newDipoleSource); 
	}

	public DataTag getTag() {
		return tag;
	}

	public IEtomicaDataInfo getDataInfo() {
		return dataInfo;
	}

	public MoleculeAgent makeAgent(IAtom a,IBox box) {
		 return new MoleculeAgent(space);
	}

	public void releaseAgent(MoleculeAgent agent, IAtom a,IBox box) {
		
	}

	public Class getAgentClass() {
		return MoleculeAgent.class;
	}

    public static class MoleculeAgent implements Integrator.Torquable, Integrator.Forcible, Serializable {  //need public so to use with instanceof
        private static final long serialVersionUID = 1L;
        public final IVectorMutable torque;
        public final IVectorMutable force;
        
        public MoleculeAgent(ISpace space) {
            torque = new Vector1D();
            force = space.makeVector();
        }
        
        public IVectorMutable torque() {return torque;}
        public IVectorMutable force() {return force;}
    }
	
}
