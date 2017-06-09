/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.paracetamol;

import etomica.box.Box;
import etomica.atom.IMoleculeList;
import etomica.potential.PotentialMaster;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.normalmode.CoordinateDefinition;
import etomica.math.function.FunctionMultiDimensionalDifferentiable;
import etomica.math.numerical.FiniteDifferenceDerivative;


public class NumericalDerivativeEnergyParacetamol implements FunctionMultiDimensionalDifferentiable{
	
	public NumericalDerivativeEnergyParacetamol(Box box, PotentialMaster potentialMaster){
		this.box = box;
		this.potentialMaster = potentialMaster;
		this.meterEnergy = new MeterPotentialEnergy(potentialMaster);
		meterEnergy.setBox(box);
		
	}
	
	
	public double f(double[] u){
		
		
		for (int cell=0; cell<coordinateDefinition.getBasisCells().length; cell++){
			IMoleculeList molecules = coordinateDefinition.getBasisCells()[cell].molecules;
			coordinateDefinition.setToU(molecules, u);
		}
		
		double energy = meterEnergy.getDataAsScalar();
		
		//System.out.println("Coordinate of first atom that overlaps: "+((IAtomPositioned)((IAtomGroup)box.getMoleculeList(species).getAtom(101)).getChildList().getAtom(2)).getPosition());
        //System.out.println("Coordinate of second atom that overlaps: "+((IAtomPositioned)((IAtomGroup)box.getMoleculeList(species).getAtom(103)).getChildList().getAtom(18)).getPosition());
        
        //System.out.println(Arrays.toString(u));
		
		if (Double.isNaN(energy)|| Double.isInfinite(energy)){
			
			throw new RuntimeException("Energy is "+ energy);
		}
		return energy;
	}
	
	public double df(int[] d, double[] u){
		
		for (int cell=0; cell<coordinateDefinition.getBasisCells().length; cell++){
			IMoleculeList molecules = coordinateDefinition.getBasisCells()[cell].molecules;
			coordinateDefinition.setToU(molecules, u);
		}
		
		finiteDifferenceDerivative = new FiniteDifferenceDerivative(this);
		finiteDifferenceDerivative.setH(0.00001);
		
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

		fPrime[index] = finiteDifferenceDerivative.df(d, u);
		return fPrime[index];
	
	}
	
	public int getDimension(){
		return coordinateDefinition.getCoordinateDim();
	}
	
	public void setCoordinateDefinition(CoordinateDefinition coordinateDefinition){
		this.coordinateDefinition = coordinateDefinition;
		fPrime = new double[coordinateDefinition.getCoordinateDim()];
	}
	
	public CoordinateDefinition getCoordinateDefinition(){
		return this.coordinateDefinition;
	}
	
	
	protected FiniteDifferenceDerivative finiteDifferenceDerivative;
	protected CoordinateDefinition coordinateDefinition;
	protected MeterPotentialEnergy meterEnergy;
	protected double[] fPrime;
	protected Box box;
	protected PotentialMaster potentialMaster;
	//public Species species;
	
}
