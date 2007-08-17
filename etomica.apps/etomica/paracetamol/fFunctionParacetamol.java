package etomica.paracetamol;

import java.io.Serializable;

import etomica.action.Activity;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomSet;
import etomica.atom.IAtomGroup;
import etomica.atom.IAtomPositioned;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntegratorHardField.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.IVector;
import etomica.space3d.IVector3D;

public class fFunctionParacetamol extends etomica.conjugategradient.fFunction implements Serializable{
	
	protected Box box;
	protected MeterPotentialEnergy meterEnergy;
	protected PotentialMaster potentialMaster;
	protected IteratorDirective allAtoms;
	protected PotentialCalculationForceSum forceSum;
	protected AtomAgentManager agentManager;
	protected Activity activity;
	protected int coordinateDim;
	
	public fFunctionParacetamol(Box box, PotentialMaster potentialMaster){
		super(box, potentialMaster);
		rotationAxis = (IVector3D)box.getSpace().makeVector();
		a      = (IVector3D)box.getSpace().makeVector();
		aProj  = (IVector3D)box.getSpace().makeVector();
		v      = (IVector3D)box.getSpace().makeVector();
		deltaV = (IVector3D)box.getSpace().makeVector();
		d      = new IVector3D[20];
		torque = new IVector3D[20];
		torqueF= new IVector3D[20];
		torqueSum = new IVector3D[8]; // only looking at one cell
	}
	
	public double f(){
		super.f();
		return meterEnergy.getDataAsScalar();
	}
	
	public double[] fPrime(AtomSet molecules, double[] u){ 
		
		/*
		 * u is the generalized coordinate
		 */
		
		
		// u = the position with delta or the specified position
		// need to call setToU with the u values
		
		coordDef.setToU(molecules,u);
		
		double[] fPrime = super.fPrime(molecules, u);
		
		forceSum.reset();
	
		/*
		 *  fPrime[coordinateDim] 
		 * 	where we have 6 generalized coordinates: 3 modes on translation and 3 on rotation for each molecule
		 *  
		 */
		
		int j=3;
		
		for (int p=0; p<molecules.getAtomCount(); p++){
			IAtomGroup molecule = (IAtomGroup)molecules.getAtom(p);
		
			 //leafPos0 is atom C1 in Paracetamol
			 //leafPos5 is atom C4 in Paracetamol
			IVector leafPos0 = ((IAtomPositioned)molecule.getChildList().getAtom(0)).getPosition();
			IVector leafPos5 = ((IAtomPositioned)molecule.getChildList().getAtom(5)).getPosition();
			
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
				 }
				 
				 if(jay==2){
					 rotationAxis.E(v);
				 }
				 
				 
				 /*
				  * To find the distance vector of each atoms within p-th molecule
				  * that is perpendicular to the rotation axis
				  */
				 for (int q=0; q<molecule.getChildList().getAtomCount(); q++){
					 
		    	    	/*
		    	    	 * Determine the distance, d, by using Vector Projection
		    	    	 */
		    	    	a.Ev1Mv2(((IAtomPositioned)molecule.getChildList().getAtom(q)).getPosition(), leafPos0);
		    	    	
		    	    	double dotProd = a.dot(rotationAxis);
		    	    	aProj.Ea1Tv1(dotProd,rotationAxis);
		    	    	
		    	    	d[q].Ev1Mv2(a, aProj);
		    	    	d[q].normalize();
				 }
				 
				 /*
				  * The forces acting on each individual atoms within a given p-th molecule
				  *   x-component, y-componet and z-component
				  *   
				  *   And summing the torque of all atoms to torqueSum[]
				  */
				 
				 for (int q=0; q<molecule.getChildList().getAtomCount(); q++){ 
					
					deltaV.E(d[q]);
					deltaV.XE(rotationAxis);
					deltaV.normalize();
				
					double scalarF = ((IntegratorVelocityVerlet.MyAgent)agentManager.getAgent(molecules.getAtom(q))).force()
										.dot(deltaV);
					torqueF[q].Ea1Tv1(scalarF, deltaV);
					torque[q].E(d[q]);                         // torque = d X F
					torque[q].XE(torqueF[q]);
					
					torqueSum[p].PE(torque[q]);    
					
				 }
				 
				fPrime[j+jay] = Math.sqrt(torqueSum[p].squared());  //taking the magnitude
				
			 }
			
			 j += coordinateDim/molecules.getAtomCount();
		}
		
		return fPrime;
	}
	
	protected final IVector3D rotationAxis;
	protected final IVector3D a, aProj;
	protected final IVector3D v, deltaV;
	protected final IVector3D [] d, torque, torqueF, torqueSum;
	public CoordinateDefinitionParacetamol coordDef;
	private static final long serialVersionUID = 1L;
}
