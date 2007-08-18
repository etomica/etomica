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
import etomica.normalmode.CoordinateDefinition;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;

public class fFunction {

	/*
	 * This class is developed to calculate for the first derivative of energy function.
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
	protected CoordinateDefinition coordinateDefinition;
	protected double[] fPrime;
	
	public fFunction(Box box, PotentialMaster potentialMaster){
		this.box = box;
		this.potentialMaster = potentialMaster;
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
	
	public CoordinateDefinition getCoordinateDefinition(){
		return coordinateDefinition;
	}
	
	public void setCoordinateDefinition (CoordinateDefinition coordinateDefinition){
		this.coordinateDefinition = coordinateDefinition;
		fPrime = new double[coordinateDefinition.getCoordinateDim()];
	}
	
	public double[] fPrime(double[] newU){
		
		/*
		 * Index is assigned to be the number of molecules in a box
		 * fPrime[number]; number equals the degree of freedom, 
		 * 	where dx, dy, and dz to each molecule
		 */
		
		forceSum.reset();
		
		int j=0;
		potentialMaster.calculate(box, allAtoms, forceSum);
		
		AtomSet molecules = coordinateDefinition.getBasisCells()[0].molecules;
		
		for (int m=0; m<molecules.getAtomCount(); m++){
			
			if (m==0){
				for (int k=0;k<3;k++){
					fPrime[j+k] = 0;
				}
				
			} else
				
			for (int k=0; k<3; k++){
				fPrime[j+k] = ((IntegratorVelocityVerlet.MyAgent)agentManager.getAgent(molecules.getAtom(m)))
							   .force.x(k);
			}
		
			j += coordinateDefinition.getCoordinateDim() /molecules.getAtomCount();
		}
		
		return fPrime;
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
