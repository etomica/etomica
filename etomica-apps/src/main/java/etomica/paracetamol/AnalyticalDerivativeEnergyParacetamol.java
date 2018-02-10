/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.paracetamol;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.conjugategradient.DerivativeEnergyFunction;
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;

import java.io.Serializable;

public class AnalyticalDerivativeEnergyParacetamol extends DerivativeEnergyFunction implements Serializable{
	
	public AnalyticalDerivativeEnergyParacetamol(Box box, PotentialMaster potentialMaster, Space space){
		super(box, potentialMaster, space);
		rotationAxis = space.makeVector();
		a      = space.makeVector();
		aProj  = space.makeVector();
		v      = space.makeVector();
		deltaV = space.makeVector();
		distance = new Vector[20];
		torque   = new Vector[20];
		torqueF  = new Vector[20];
		for (int i=0; i<20; i++){
			distance[i] = space.makeVector();
			torque  [i] = space.makeVector();
			torqueF [i] = space.makeVector();
		}
		torqueSum = space.makeVector();
		
	}
	
	public double df(int[] d, double[] u){ 
		
		/*
		 * u is the generalized coordinate
		 */
		
		/*
		 * d only takes in array that compute first-order derivative w.r.t. to corresponding n-th dimension
		 *  for example, d=new double{1, 0, 0} or {0, 0, 1}, which means first-order differentiation to 
		 *  first- and third- dimension respectively. 
		 */
		
		int index =0;
		double check =0;
		fPrimeRotation = new double[u.length];
		
		for (int i =0; i <d.length; i++){
			check += d[i];
			
			if (d[i]==1){
				index = i;
			}
		} 
		
		if (check != 1){
			throw new IllegalArgumentException("The function MUST and CAN only compute first-order derivative!!");
		}
		
		int group = index;
		
		/*
		 * return the first-derivative of translation mode
		 */
		if(group%6 < 3){
			
			fPrimeRotation[index] = super.df(d, u);
			return fPrimeRotation[index];
			
		} else {
			
			forceSum.reset();

			for (int cell=0; cell<coordinateDefinition.getBasisCells().length; cell++){
				IMoleculeList molecules = coordinateDefinition.getBasisCells()[cell].molecules;
				coordinateDefinition.setToU(molecules, u);
			}
			
			/*
			 *  fPrime[coordinateDim] 
			 * 	where we have 6 generalized coordinates: 3 modes on translation and 3 on rotation for each molecule
			 *  
			 */
			IMoleculeList molecules = coordinateDefinition.getBasisCells()[0].molecules;
			
			int j=3;
			
			for (int p=0; p<molecules.getMoleculeCount(); p++){ //loop over the 8 molecules in the basis cell
				
				IAtomList molecule = molecules.getMolecule(p).getChildList();
			
				 //leafPos0 is atom C1 in Paracetamol
				 //leafPos5 is atom C4 in Paracetamol
				Vector leafPos0 = molecule.getAtom(0).getPosition();
				Vector leafPos5 = molecule.getAtom(5).getPosition();
				
				v.Ev1Mv2(leafPos5, leafPos0);
				v.normalize();
				 
				potentialMaster.calculate(box, allAtoms, forceSum);
				 
				/*
				 * To find the rotation axis by taking the cross-product
				 * of v and delta v
				 * 
				 * there are 3 cases: u[3], u[4], and u[5]
				 * setting the rotation axis correspondingly
				 * 
				 *  having jay to loop over the 3 cases
				 */
				 
				for (int jay=0; jay<3; jay++){
					 
					 if(jay==0){
					 	deltaV.E(new double[]{-u[j]/Math.sqrt(1-u[j]*u[j]-u[j+1]*u[j+1]) ,1 ,0});
					 	rotationAxis.E(v);
					 	rotationAxis.XE(deltaV);
					 	rotationAxis.normalize();
					 } else
					
					 if(jay==1){
						 deltaV.E(new double[]{-u[j+1]/Math.sqrt(1-u[j]*u[j]-u[j+1]*u[j+1]) ,0 ,1});
						 rotationAxis.E(v);
						 rotationAxis.XE(deltaV);
						 rotationAxis.normalize();
					 } else
					 
					 if(jay==2){
						 rotationAxis.E(v);
					 }
					 
					 
					 /*
					  * To find the distance vector, d[] of each atoms within p-th molecule
					  * that is perpendicular to the rotation axis
					  */
					 for (int q=0; q<molecule.getAtomCount(); q++){
						 
			    	    	/*
			    	    	 * Determine the distance, d, by using Vector Projection
			    	    	 */
						 
						 	// vector a when q=0
						 	if (q==0){
						 		a.E(new double[] {0, 0, 0});
						 	} else {
						 	
				    	    	a.Ev1Mv2(molecule.getAtom(q).getPosition(), leafPos0);
				    	    	a.normalize();
						 	}
						 	
			    	    	double dotProd = a.dot(rotationAxis);
			    	    	aProj.Ea1Tv1(dotProd,rotationAxis);
			    	    	aProj.normalize();
			    	    	
			    	    	if (q==0){
						 		distance[q].E(new double[] {0, 0, 0});
						 	} else {
						 		
						 		distance[q].Ev1Mv2(a, aProj);             
						 		distance[q].normalize();                 
						 	}
					 }
					 
					 /*
					  * The forces acting on each individual atoms within a given p-th molecule
					  *   x-component, y-componet and z-component
					  *   
					  *   And summing the torque of all atoms to torqueSum[]
					  */
	
					 moleculeForce.E(0); //initialize moleculeForce to zero
					 torqueSum.E(0);
					 
					 for (int q=0; q<molecule.getAtomCount(); q++){ 
						
						if (q==0){
							deltaV.E(new double[] {0, 0, 0});
						} else {
							
						deltaV.E(distance[q]);
						deltaV.XE(rotationAxis);
						deltaV.normalize();
						}
						
						moleculeForce.E(agentManager.getAgent(molecule.getAtom(q)));
						
						double scalarF = 0;
						scalarF = moleculeForce.dot(deltaV);
						torqueF[q].Ea1Tv1(scalarF, deltaV);
						
						if (q==0){
					 		torque[q].E(new double[] {0, 0, 0});
					 	} else {
							torque[q].E(d[q]);                         // torque = d X F
							torque[q].XE(torqueF[q]);
					 	}
						
						torqueSum.PE(torque[q]);    
						// torqueSum all equal to NaN!!!!!!!!!!!!!
					 }
					 
					fPrimeRotation[j+jay] = Math.sqrt(torqueSum.squared());  //taking the magnitude
					
					if (index == j+jay){
						return fPrimeRotation[j+jay];
					}
					
				 }
				
				 j += coordinateDefinition.getCoordinateDim()/molecules.getMoleculeCount();
			}
		
		return fPrimeRotation[index];
		}
	}
	
	
	
	protected final Vector rotationAxis;
	protected final Vector a, aProj;
	protected final Vector v, deltaV;
	protected final Vector[] distance, torque, torqueF;
	protected final Vector torqueSum;
	protected double[] fPrimeRotation;
	private static final long serialVersionUID = 1L;
}
