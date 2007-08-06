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
	 * derivative of energy function.
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
	
	public fFunction(Box box, PotentialMaster potentialMaster){
		this.box = box;
		this.potentialMaster = potentialMaster;
		this.coordinateDim = 48;
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
		
		int j=0;
		potentialMaster.calculate(box, allAtoms, forceSum);
		
		for (int m=0; m<molecules.getAtomCount(); m++){
			if (m==0){
				for (int k=0;k<3;k++){
					fPrime[j+k] = 0;
				}
			} else {
				
				for (int k=0; k<3; k++){
					fPrime[j+k] =((IntegratorVelocityVerlet.MyAgent)agentManager.getAgent(molecules.getAtom(m)))
							.force.x(k);
				}
			}
			j+=coordinateDim /molecules.getAtomCount();
		}
		return fPrime;
	}
	
	/*
	 * There are 6 modes; with 8 molecules within a cell
	 * for 2-D array, it will be a 48 X 48 matrix
	 * 
	 */

	
	public double[][] fDoublePrime(AtomSet molecules, double[] newU){
		double[][] fDoublePrime = new double [coordinateDim][coordinateDim]; 
		
		//loop over the molecules within the cell
		
		for (int p=0; p<molecules.getAtomCount(); p++){
		
			for(int q=0; q<coordinateDim; q++){
				
				fDoublePrime[p][q] = fPrime(molecules)[q];
			}
		}
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
