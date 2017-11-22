/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.conjugategradient;

import etomica.action.Activity;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.math.function.FunctionMultiDimensionalDifferentiable;
import etomica.molecule.IMoleculeList;
import etomica.normalmode.CoordinateDefinition;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;

public class DerivativeEnergyFunction implements FunctionMultiDimensionalDifferentiable{

	/*
	 * This class is developed to calculate for 
	 * 	the first derivative of energy function ANALYTICALLY.
	 * 
	 * @author Tai Tan
	 */
	
	protected Box box;
	protected MeterPotentialEnergy meterEnergy;
	protected PotentialMaster potentialMaster;
	protected IteratorDirective allAtoms;
	protected PotentialCalculationForceSum forceSum;
	protected AtomLeafAgentManager<IntegratorVelocityVerlet.MyAgent> agentManager;
	protected Activity activity;
	protected CoordinateDefinition coordinateDefinition;
	protected double[] fPrime;
	protected Vector moleculeForce;
	protected FunctionMultiDimensionalDifferentiable fFunction;
	
	public DerivativeEnergyFunction(Box box, PotentialMaster potentialMaster, Space space){
		this.box = box;
		this.potentialMaster = potentialMaster;
		meterEnergy = new MeterPotentialEnergy(potentialMaster);
		allAtoms = new IteratorDirective();
		forceSum = new PotentialCalculationForceSum();
		
		MyAgentSource source = new MyAgentSource(space);
		agentManager = new AtomLeafAgentManager<IntegratorVelocityVerlet.MyAgent>(source, box);
		forceSum.setAgentManager(agentManager);
		moleculeForce = space.makeVector();
		
		/*
		 * Dimensions of a study system is three-times the number of atoms
		 */
	}
	

	
	public CoordinateDefinition getCoordinateDefinition(){
		return coordinateDefinition;
	}
	
	public void setCoordinateDefinition (CoordinateDefinition coordinateDefinition){
		this.coordinateDefinition = coordinateDefinition;
		fPrime = new double[coordinateDefinition.getCoordinateDim()];
	}
	
	public double f(double[] newU){
		for (int cell=0; cell<coordinateDefinition.getBasisCells().length; cell++){
			IMoleculeList molecules = coordinateDefinition.getBasisCells()[cell].molecules;
			coordinateDefinition.setToU(molecules, newU);
		}
		
		return meterEnergy.getDataAsScalar();
	}
	
	public double df(int[] d, double[] u){
		
		/*
		 * Index is assigned to be the number of molecules in a box
		 * fPrime[number]; number equals the degree of freedom, 
		 * 	where dx, dy, and dz to each molecule
		 */
		
		forceSum.reset();
		
		for (int cell=0; cell<coordinateDefinition.getBasisCells().length; cell++){
			IMoleculeList molecules = coordinateDefinition.getBasisCells()[cell].molecules;
			coordinateDefinition.setToU(molecules, u);
		}
		
		/*
		 * d only takes in array that compute first-order derivative w.r.t. to corresponding n-th dimension
		 *  for example, d=new double{1, 0, 0} or {0, 0, 1}, which means first-order differentiation to 
		 *  first- and third- dimension respectively. 
		 */
		
		int index =0;
		double check =0;
		
		for (int i =0; i <d.length; i++){
			check += d[i];
			
			if (d[i]==1){
				index = i;
			}
		} 
		
		if (check != 1){
			throw new IllegalArgumentException("The function MUST and CAN only compute first-order derivative!!");
		}
		
		int j=0;
		potentialMaster.calculate(box, allAtoms, forceSum);
		
		IMoleculeList molecules = coordinateDefinition.getBasisCells()[0].molecules;
		
		for (int m=0; m<molecules.getMoleculeCount(); m++){
				
			if (m==0){
				for (int k=0; k<3; k++){
					fPrime[j+k] = 0;
					
					if (index == j+k){
						return fPrime[j+k];
					}
				}
					
			} else {
				
				IAtomList childList = molecules.getMolecule(m).getChildList();
				
				moleculeForce.E(0); //initialize moleculeForce to zero
				
				for (int r=0; r<childList.getAtomCount(); r++){
					moleculeForce.PE(agentManager.getAgent(childList.getAtom(r)).force);
				}
					
				
				for (int k=0; k<3; k++){
					fPrime[j+k] = moleculeForce.getX(k);
						
					if (index == j+k){
						return fPrime[j+k];
					}
				}
			}
			j += coordinateDefinition.getCoordinateDim() /molecules.getMoleculeCount();
		}
		
		return fPrime[index];
	}
	
	public int getDimension(){
		return fFunction.getDimension();
	}
	
	public void getScalarEnergy(){
		meterEnergy.setBox(box);
		System.out.println("The energy of the system is: "+meterEnergy.getDataAsScalar());
	}
	
	
	
	public static class MyAgentSource implements AgentSource<IntegratorVelocityVerlet.MyAgent> {
		
		public MyAgentSource(Space space){
			this.space = space;
		}
		
		public void releaseAgent(IntegratorVelocityVerlet.MyAgent agent, IAtom atom, Box agentBox){}

		public IntegratorVelocityVerlet.MyAgent makeAgent(IAtom atom, Box agentBox){
			
		    return new IntegratorVelocityVerlet.MyAgent(space);
		}
		protected Space space;
	}

}
