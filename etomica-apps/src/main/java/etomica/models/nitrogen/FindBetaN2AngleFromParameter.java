/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.atom.IMolecule;
import etomica.space.Vector;
import etomica.box.Box;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.normalmode.BasisBigCell;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Degree;



/**
 * Class that creates:
 *  1. alpha angles
 *  2.  beta angles 
 *  3. rotation axis
 *  
 *  for CoordinateDefinitionNitrogen setOrientationVectorBetaLatticeSum Method
 * 
 * @author Tai Boon Tan
 *
 */
public class FindBetaN2AngleFromParameter extends Simulation{

	public FindBetaN2AngleFromParameter(Space space, double density, double[][]param) {
		super(space);
		double ratio = 1.631;
		double a = Math.pow(4.0/(Math.sqrt(3.0)*ratio*density), 1.0/3.0);
		double c = a*ratio;
		
		int nC = 2;
		
		Basis basisHCP = new BasisHcp();
		Basis basis = new BasisBigCell(space, basisHCP, new int[]{nC, nC, nC});
		
		SpeciesN2 species = new SpeciesN2(space);
		addSpecies(species);
		
		Box box = new Box(space);
		addBox(box);
		box.setNMolecules(species, nC*nC*nC*2);		
		int [] nCells = new int[]{1,1,1};
		
		Vector[] boxDim = new Vector[3];
		boxDim[0] = space.makeVector(new double[]{nC*a, 0, 0});
		boxDim[1] = space.makeVector(new double[]{-nC*a*Math.cos(Degree.UNIT.toSim(60)), nC*a*Math.sin(Degree.UNIT.toSim(60)), 0});
		boxDim[2] = space.makeVector(new double[]{0, 0, nC*c});
		
		Primitive primitive = new PrimitiveHexagonal(space, (nC)*a, nC*c);
		
		CoordinateDefinitionNitrogen coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordinateDef.setIsBeta();
		coordinateDef.setOrientationVectorBeta(space);
		coordinateDef.initializeCoordinates(nCells);

		double[] u = new double[20];
		
		int kParam=0;
		for (int i=0; i<param.length;i++){
			for (int j=0; j<param[0].length;j++){
				u[kParam]=param[i][j];
				kParam++;
			}	
		}
					
		int numDOF = coordinateDef.getCoordinateDim();
		double[] newU = new double[numDOF];
		
		for(int j=0; j<numDOF; j+=10){
			if(j>0 && j%(nC*10)==0){
				j+=nC*10;
				if(j>=numDOF){
					break;
				}
			}
			for(int k=0; k<10;k++){
				newU[j+k]= u[k];
			}
		}
		
		for(int j=nC*10; j<numDOF; j+=10){
			if(j>nC*10 && j%(nC*10)==0){
				j+=nC*10;
				if(j>=numDOF){
					break;
				}
			}
			for(int k=0; k<10;k++){
				newU[j+k]= u[k+10];
			}
		}
		
		coordinateDef.setToU(box.getMoleculeList(), newU);
		coordinateDef.initNominalU(box.getMoleculeList());
		
		Vector[] aVector = new Vector[4];
		Vector[] cVector = new Vector[4];
		rotationAxis = new Vector[4];
		deviationVector = new Vector[4];
		
		Vector temp1 = space.makeVector();
		Vector temp2 = space.makeVector();
		Vector bVector = space.makeVector(new double[]{0.0, 0.0, 1.0});
		
		for (int i=0; i<aVector.length; i++){
			aVector[i] = space.makeVector();
			cVector[i] = space.makeVector();
			rotationAxis[i] = space.makeVector();
			deviationVector[i] = space.makeVector();
		}
		
		alpha = new double[4];
		beta = new double[4];
		
		for(int i=0; i<4; i++){
			int j;
			if(i<2){
				j=i;
			} else {
				j=i+2;
			}
			
			IMolecule molecule = coordinateDef.getBox().getMoleculeList().getMolecule(j);
		  	Vector molleafPos0 = molecule.getChildList().getAtom(0).getPosition();
		   	Vector molleafPos1 = molecule.getChildList().getAtom(1).getPosition();
		   	
			aVector[i].Ev1Mv2(molleafPos1, molleafPos0);
		    aVector[i].normalize();
			
		    temp1.E(aVector[i]);
		    temp1.XE(bVector);
		    temp1.normalize();
		    temp2.E(bVector);
		    temp2.XE(temp1);
		    cVector[i].Ea1Tv1(Math.sqrt(aVector[i].squared()), temp2);
		    cVector[i].normalize();
		    
		    temp1.E(new double[]{1.0, 0.0, 0.0});
		    alpha[i] = Math.acos(cVector[i].dot(temp1));
		    beta[i]  = Math.acos(aVector[i].dot(cVector[i]));
		    
		}
		
		int contA, contB;
//		System.out.println("rho: " + density);
	    for(int i=0; i<aVector.length; i++){
	    	contA = contB = 1;
	    	
	    	if(i==0){contB=-1;}
	    	else if (i==2){contA=-1;}
	    	else if(i==3){contA=-1;contB=-1;}
	    	
	    	alpha[i] = contA*Degree.UNIT.fromSim(alpha[i]);
	    	beta [i] = contB*Degree.UNIT.fromSim(beta[i]);
	    	
	    	rotationAxis[i].E(new double[]{ cVector[i].getX(1), -cVector[i].getX(0), 0.0});
		    deviationVector[i].E(new double[]{param[i][0], param[i][1], param[i][2]});
	    	
//	    	System.out.println("aVector["+i+"]   : " + aVector[i].toString());
//	    	System.out.println("cVector["+i+"]   : " + cVector[i].toString());
//	    	System.out.println("alpha["+i+"](deg): " + alpha[i]);
//	    	System.out.println(" beta["+i+"](deg): " + beta[i]);
//	    	System.out.println("rotation["+i+"]: " + rotationAxis[i]);
//		    System.out.println("xyz:["+i+"]: "+ deviationVector[i].toString());
//	    	System.out.println();
	    }
	 
	}
	
	public Vector[] getDeviationVector() {
		return deviationVector;
	}

	public double[] getAlpha() {
		return alpha;
	}

	public double[] getBeta() {
		return beta;
	}

	public Vector[] getRotationAxis() {
		return rotationAxis;
	}

	public static void main (String[] args){
		
		int[] nC = new int []{2,2,2};
		double density = 0.0230;
		BetaPhaseLatticeParameter parameters = new BetaPhaseLatticeParameter();
		double[][] param = parameters.getParameter(density);
		
		final FindBetaN2AngleFromParameter sim = new FindBetaN2AngleFromParameter(Space3D.getInstance(3), density, param);

	}

	protected double[] alpha, beta;
	protected Vector[] rotationAxis;
	protected Vector[] deviationVector;
	private static final long serialVersionUID = 1L;
}
