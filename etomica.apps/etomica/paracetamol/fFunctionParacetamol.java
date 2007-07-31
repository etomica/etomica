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
	
	public fFunctionParacetamol(Box box, PotentialMaster potentialMaster, int coordinateDim){
		super(box, potentialMaster, coordinateDim);
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
		
		return meterEnergy.getDataAsScalar();
	}
	
	public double[] fPrime(){
		
		forceSum.reset();
		/*
		 * Index is assigned to be the number of molecules in a box
		 * fPrime[number]; number equals the degree of freedom, 
		 * 	where dx, dy, and dz to each molecule
		 */
		int index = box.getLeafList().getAtomCount();
		double[] fPrime = new double [coordinateDim];
		
		int j=3;
		potentialMaster.calculate(box, allAtoms, forceSum);
		
		
		for (int m=0; m<index; m++){
			for (int k=0; k<3; k++){
				fPrime[j] =((IntegratorVelocityVerlet.MyAgent)agentManager.getAgent(box.getLeafList().getAtom(m)))
						.force.x(k);
				j+=coordinateDim/box.getLeafList().getAtomCount();
			}
		}
		return fPrime;
	}
	
	public double[] fPrime(AtomSet molecules, double[] newU){
		
		super.fPrime(molecules);
		
		forceSum.reset();
		/*
		 * Index is assigned to be the number of molecules in a box
		 * fPrimeMode[number]; number equals 6-times Index, 
		 * 	where we have 6 modes: 3 modes on translation and 3 on rotation for each molecule
		 */
		
		
		double[] fPrime = new double [coordinateDim];
		
		int j=3;
		
		/*
		 * Picking vectors to define the rotational angle(theta)
		 */
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
			  * Creating 3 cases: u[3], u[4], and u[5]
			  * setting the rotation axis correspondingly
			  */
			 
			 //for loop for looping over j
			 
			 for (int jay=j; jay<j+3; jay++){
				 
				 if(jay%6 ==3){
				 	deltaV.assignTo(new double[]{-newU[j]/Math.sqrt(1-newU[j]*newU[j]-newU[j+1]*newU[j+1]) ,1 ,0});
				 	rotationAxis.E(v);
				 	rotationAxis.XE(deltaV);
				 	rotationAxis.normalize();
				 } else
				
				 if(jay%6 ==4){
					 deltaV.assignTo(new double[]{-newU[j]/Math.sqrt(1-newU[j-1]*newU[j-1]-newU[j]*newU[j]) ,0 ,1});
					 rotationAxis.E(v);
					 rotationAxis.XE(deltaV);
					 rotationAxis.normalize();
				 }
				 
				 if(jay%6 ==5){
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
		    	    	aProj.Ea1Tv1(dotProd/rotationAxis.squared(),rotationAxis);
		    	    	
		    	    	d[q].Ev1Mv2(a, aProj);
		    	    	d[q].normalize();
		    	    	
				 }
				 
				 /*
				  * The forces acting on each invidual atoms within a given p-th molecule
				  *   x-component, y-componet and z-component
				  *   
				  *   And summing the torque of all atoms to torqueSum[]
				  */
				 for (int q=0; q<molecule.getChildList().getAtomCount(); q++){ 
					 
					IVector3D []force = new IVector3D [molecule.getChildList().getAtomCount()];
				
					if (jay%6==5){
						deltaV.E(d[10]);    //distance of atom10 from rotation axis
						deltaV.XE(rotationAxis);
						deltaV.normalize();
					}
					
					force[q].E(((IntegratorVelocityVerlet.MyAgent)agentManager.getAgent(molecules.getAtom(q))).force());
					double scalarF = force[q].dot(deltaV);
					torqueF[q].Ea1Tv1(scalarF, deltaV);
					torque[q].E(d[q]);                         // torque = r X F
					torque[q].XE(torqueF[q]);
					
					torqueSum[p].PE(torque[q]);              
					
				 }
				 
				 
				 if (jay%6==5){
					 fPrime[jay] = torqueSum[p].x(0)+torqueSum[p].x(1)+torqueSum[p].x(2); ///********************check torqueSum[p].dot(1,1,1)
				 } else
					 fPrime[jay]  = torqueSum[p].dot(deltaV);
			 }
			 
			 j += coordinateDim/molecules.getAtomCount();
		}
		
		return fPrime;
	}
	
	public double[][] fDoublePrime(){
		int index = box.getLeafList().getAtomCount();
		double[][] fDoublePrime = new double [9*index][9*index]; 
		
		
		
		return fDoublePrime;
	}
	
	
	protected final IVector3D rotationAxis;
	protected final IVector3D a, aProj;
	protected final IVector3D v, deltaV;
	protected final IVector3D [] d, torque, torqueF, torqueSum;
	private static final long serialVersionUID = 1L;
}
