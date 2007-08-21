package etomica.conjugategradient;

import etomica.action.Activity;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.IAtomGroup;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.normalmode.CoordinateDefinition;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.util.FunctionMultiDimensionalDifferentiable;

public class DerivativeFunction implements FunctionMultiDimensionalDifferentiable {

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
	protected IVector moleculeForce;
	
	public DerivativeFunction(Box box, PotentialMaster potentialMaster){
		this.box = box;
		this.potentialMaster = potentialMaster;
		meterEnergy = new MeterPotentialEnergy(potentialMaster);
		allAtoms = new IteratorDirective();
		forceSum = new PotentialCalculationForceSum();
		
		MyAgentSource source = new MyAgentSource(box.getSpace());
		agentManager = new AtomAgentManager(source, box);
		forceSum.setAgentManager(agentManager);
		moleculeForce = box.getSpace().makeVector();
		
		/*
		 * Dimensions of a study system is three-times the number of atoms
		 */
	}
	
	public double function(double[] newU){
		for (int cell=0; cell<coordinateDefinition.getBasisCells().length; cell++){
			AtomSet molecules = coordinateDefinition.getBasisCells()[cell].molecules;
			coordinateDefinition.setToU(molecules, newU);
		}
		
		return meterEnergy.getDataAsScalar();
	}
	
	public CoordinateDefinition getCoordinateDefinition(){
		return coordinateDefinition;
	}
	
	public void setCoordinateDefinition (CoordinateDefinition coordinateDefinition){
		this.coordinateDefinition = coordinateDefinition;
		fPrime = new double[coordinateDefinition.getCoordinateDim()];
	}
	
	public double[] dfdx(double[] newU){
		
		/*
		 * Index is assigned to be the number of molecules in a box
		 * fPrime[number]; number equals the degree of freedom, 
		 * 	where dx, dy, and dz to each molecule
		 */
		
		forceSum.reset();
		
		for (int cell=0; cell<coordinateDefinition.getBasisCells().length; cell++){
			AtomSet molecules = coordinateDefinition.getBasisCells()[cell].molecules;
			coordinateDefinition.setToU(molecules, newU);
		}
		
		int j=0;
		potentialMaster.calculate(box, allAtoms, forceSum);
		
		AtomSet molecules = coordinateDefinition.getBasisCells()[0].molecules;
		
		for (int m=0; m<molecules.getAtomCount(); m++){
			
			if (m==0){
				for (int k=0;k<3;k++){
					fPrime[j+k] = 0;
				}
				
			} else {
				
				AtomSet childList = ((IAtomGroup)molecules.getAtom(m)).getChildList();
				
				moleculeForce.E(0); //initialize moleculeForce to zero
				
				for (int r=0; r<childList.getAtomCount(); r++){
					moleculeForce.PE(((IntegratorVelocityVerlet.MyAgent)agentManager.getAgent(childList.getAtom(r)))
							   .force);
				}
				
			
				for (int k=0; k<3; k++){
					fPrime[j+k] = moleculeForce.x(k);
				}
			}
			j += coordinateDefinition.getCoordinateDim() /molecules.getAtomCount();
		}
		
		return fPrime;
	}
	
	public void getScalarEnergy(){
		meterEnergy.setBox(box);
		System.out.println("The energy of the system is: "+meterEnergy.getDataAsScalar());
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
