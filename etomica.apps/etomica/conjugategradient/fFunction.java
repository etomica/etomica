package etomica.conjugategradient;

import etomica.action.Activity;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntegratorHardField.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;

public class fFunction {

	/*
	 * This class is developed to calculate for the first and second
	 * derivation of energy function.
	 * 
	 * @author Tai Tan
	 */
	
	protected Box box;
	protected MeterPotentialEnergy meterEnergy;
	protected PotentialMaster potentialMaster;
	protected IteratorDirective allAtoms;
	protected PotentialCalculationForceSum forceSum;
	protected AtomAgentManager agentManager;
	protected Activity activity;
	protected int coordinateDim;
	
	public fFunction(Box box, PotentialMaster potentialMaster, int coordinateDim){
		this.box = box;
		this.potentialMaster = potentialMaster;
		this.coordinateDim = coordinateDim;
		meterEnergy = new MeterPotentialEnergy(potentialMaster);
		allAtoms = new IteratorDirective();
		forceSum = new PotentialCalculationForceSum();
		
		
		MyAgentSource source = new MyAgentSource(box.getSpace());
		agentManager = new AtomAgentManager(source, box);
		forceSum.setAgentManager(agentManager);
		/*
		 * Dimensions of a study system is three-times the number of atoms
		 */
	}
	
	public double f(){
		return meterEnergy.getDataAsScalar();
	}
	
	
	public double[] fPrime(AtomSet molecules){
		
		forceSum.reset();
		/*
		 * Index is assigned to be the number of molecules in a box
		 * fPrime[number]; number equals the degree of freedom, 
		 * 	where dx, dy, and dz to each molecule
		 */
		
		double[] fPrime = new double [coordinateDim];
		
		int j=6;
		potentialMaster.calculate(box, allAtoms, forceSum);
		
		for (int m=0; m<molecules.getAtomCount(); m++){
			for (int k=0; k<3; k++){
				fPrime[j] =((IntegratorVelocityVerlet.MyAgent)agentManager.getAgent(molecules.getAtom(m)))
						.force.x(k);
				j+=coordinateDim /molecules.getAtomCount();
			}
		}
		return fPrime;
	}
	
	public double[] fPrimeMode(AtomSet molecules, double[] newU){
		
		forceSum.reset();
		/*
		 * Index is assigned to be the number of molecules in a box
		 * fPrimeMode[number]; number equals 6-times Index, 
		 * 	where we have 6 modes: 3 modes on translation and 3 on rotation for each molecule
		 */
		int index = box.getLeafList().getAtomCount();
		double[] fPrimeMode = new double [6*index];
		
		int j=0;
		potentialMaster.calculate(box, allAtoms, forceSum);
		
		for (int p=0; p<molecules.getAtomCount(); p++){
			IAtom molecule = molecules.getAtom(p);
			
		}
		
		
		
		
		
		for (int i=0; i<index; i++){
			for (int k=0; k<3; k++){
				fPrimeMode[j] =((IntegratorVelocityVerlet.MyAgent)agentManager.getAgent(box.getLeafList().getAtom(i)))
						.force.x(k);
				j++;
			}
		}
		return fPrimeMode;
	}
	
	public double[][] fDoublePrime(){
		int index = box.getLeafList().getAtomCount();
		double[][] fDoublePrime = new double [9*index][9*index]; 
		
		
		
		return fDoublePrime;
	}
	
	public static class MyAgentSource implements AgentSource{
		
		public MyAgentSource(Space space){
			this.space = space;
		}
		
		public void releaseAgent(Object agent, IAtom atom){}
		public Class getAgentClass(){
			return IntegratorVelocityVerlet.MyAgent.class;
		}
		public Object makeAgent(IAtom atom){
			
			return new IntegratorVelocityVerlet.MyAgent(space);
			}
		protected Space space;
	}

}
