package etomica.paracetamol;

import java.io.Serializable;

import etomica.atom.AtomSet;
import etomica.atom.IAtomGroup;
import etomica.atom.IAtomPositioned;
import etomica.box.Box;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.PotentialMaster;
import etomica.space.IVector;
import etomica.space3d.IVector3D;

public class DerivativeFunctionParacetamol extends etomica.conjugategradient.DerivativeFunction implements Serializable{
	
	public DerivativeFunctionParacetamol(Box box, PotentialMaster potentialMaster){
		super(box, potentialMaster);
		rotationAxis = (IVector3D)box.getSpace().makeVector();
		a      = (IVector3D)box.getSpace().makeVector();
		aProj  = (IVector3D)box.getSpace().makeVector();
		v      = (IVector3D)box.getSpace().makeVector();
		deltaV = (IVector3D)box.getSpace().makeVector();
		d      = new IVector3D[20];
		torque = new IVector3D[20];
		torqueF= new IVector3D[20];
		for (int i=0; i<20; i++){
			d[i]      = (IVector3D)box.getSpace().makeVector();
			torque[i] = (IVector3D)box.getSpace().makeVector();
			torqueF[i]= (IVector3D)box.getSpace().makeVector();
		}
		torqueSum = box.getSpace().makeVector();
		
	}
	
	public double[] dfdx(double[] u){ 
		
		/*
		 * u is the generalized coordinate
		 */
		
		double[] fPrime = super.dfdx(u);
		
		forceSum.reset();
	
		/*
		 *  fPrime[coordinateDim] 
		 * 	where we have 6 generalized coordinates: 3 modes on translation and 3 on rotation for each molecule
		 *  
		 */
		AtomSet molecules = coordinateDefinition.getBasisCells()[0].molecules;
		
		int j=3;
		
		for (int p=0; p<molecules.getAtomCount(); p++){ //loop over the 8 atoms in the basis cell
			AtomSet molecule = ((IAtomGroup)molecules.getAtom(p)).getChildList();
		
			 //leafPos0 is atom C1 in Paracetamol
			 //leafPos5 is atom C4 in Paracetamol
			IVector leafPos0 = ((IAtomPositioned)molecule.getAtom(0)).getPosition();
			IVector leafPos5 = ((IAtomPositioned)molecule.getAtom(5)).getPosition();
			
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
					 deltaV.E(new double[]{-u[j]/Math.sqrt(1-u[j-1]*u[j-1]-u[j]*u[j]) ,0 ,1});
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
					 	
			    	    	a.Ev1Mv2(((IAtomPositioned)molecule.getAtom(q)).getPosition(), leafPos0);
			    	    	a.normalize();
					 	}
					 	
		    	    	double dotProd = a.dot(rotationAxis);
		    	    	aProj.Ea1Tv1(dotProd,rotationAxis);
		    	    	
		    	    	if (q==0){
					 		d[q].E(new double[] {0, 0, 0});
					 	} else {
					 		
					 		d[q].Ev1Mv2(a, aProj);
					 		d[q].normalize();
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
						
					deltaV.E(d[q]);
					deltaV.XE(rotationAxis);
					deltaV.normalize();
					}
					
					moleculeForce.E(((IntegratorVelocityVerlet.MyAgent)agentManager.getAgent(molecule.getAtom(q)))
								   .force);
					
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
				 
				fPrime[j+jay] = Math.sqrt(torqueSum.squared());  //taking the magnitude
				
			 }
			
			 j += coordinateDefinition.getCoordinateDim()/molecules.getAtomCount();
		}
		
		return fPrime;
	}
	
	
	
	protected final IVector3D rotationAxis;
	protected final IVector3D a, aProj;
	protected final IVector3D v, deltaV;
	protected final IVector3D [] d, torque, torqueF;
	protected final IVector torqueSum;
	private static final long serialVersionUID = 1L;
}
