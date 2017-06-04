/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.api.IMolecule;
import etomica.atom.IMoleculePositionDefinition;
import etomica.atom.MoleculePositionGeometricCenter;
import etomica.space.Vector;
import etomica.space.Space;


/**
 *  BETA-NITROGEN MODEL
 * 
 * Class to identify the pair of molecules that are the same in terms of displacement vector and orientation
 * 
 * This is done to speed up the construction of the Hessian Matrix such that we do not have to recalculate
 * 	the block matrix that is the same.
 * 
 * The index field is a 5-dimensional array
 * 	index[0] is x-direction displacement
 *  index[1] is y-direction displacement
 *  index[2] is z-direction displacement
 *  index[3] is the orientation of moleculeA 
 *  index[4] is the orientation of moleculeB
 *  
 * The isNewPair array is to keep track whether the pair with index[][][][][] 
 * has been looped over or not.
 *   
 * @author Tai Boon Tan
 *
 */
public class FindPairMoleculeIndexBetaN2 {
	public FindPairMoleculeIndexBetaN2(Space space, CoordinateDefinitionNitrogen coordinateDefination){
		this.coordinateDef = coordinateDefination;
		positionDefinition = new MoleculePositionGeometricCenter(space);
		tempVec = space.makeVector();
		tempOrientA = space.makeVector();
		tempOrientB = space.makeVector();
		
		molAVec = space.makeVector();
		molBVec = space.makeVector();
		
		int numMolecule = coordinateDefination.getBox().getMoleculeList().getMoleculeCount();
		nCell = (int) Math.round(Math.pow((numMolecule/1.999999999), 1.0/3.0));
		if(nCell > 20){
			throw new RuntimeException("<FindPairMoleculeIndexBetaNitrogen> nCell is greater than 20!!! " +
					"YOU ARE CRASHING JAVA MEMORY!! Live long and prosper!");
		}
		
		/*
		 * offset from the HCP lattice site
		 */
		latticeOffset = coordinateDefination.positionVector;
		
		//get shortest displacement vector
		molAVec.E(positionDefinition.position(coordinateDefination.getBox().getMoleculeList().getMolecule(0)));
		molBVec.E(positionDefinition.position(coordinateDefination.getBox().getMoleculeList().getMolecule(1)));
		
		molAVec.ME(latticeOffset[0]);
		molBVec.ME(latticeOffset[1]);
		
		tempVec.Ev1Mv2(molBVec, molAVec);
		
		lengthX = tempVec.getX(0);
		lengthY = tempVec.getX(1);
		lengthZ = tempVec.getX(2);
		
		int xSites = 2*nCell+1;
		int ySites = 4*nCell+1;
		int zSites = 2*nCell+1;
		siteDisplacement = new double[xSites][ySites][zSites];
		orientation = new Vector[4];
		
		for (int i=0; i<orientation.length; i++){
			orientation[i] = space.makeVector();
		}

		index = new int[5];
		isNewPair = new boolean[xSites][ySites][zSites][4][4];
		
		for(int a=0; a<isNewPair.length; a++){
			for(int b=0; b<isNewPair[0].length; b++){
				for(int c=0; c<isNewPair[0][0].length; c++){
					for(int d=0; d<isNewPair[0][0][0].length; d++){
						for(int e=0; e<isNewPair[0][0][0][0].length; e++){
							isNewPair[a][b][c][d][e] = true;
						}	
					}	
				}	
			}	
		}
	
	
	}
	
	public int[] getPairMoleculesIndex(IMolecule moleculeA, IMolecule moleculeB, boolean isReverseOrder){
	
		// Molecule A
		molAVec.E(positionDefinition.position(moleculeA));
			
		Vector molAleafPos0 = moleculeA.getChildList().getAtom(0).getPosition();
    	Vector molAleafPos1 = moleculeA.getChildList().getAtom(1).getPosition();
    	 
		tempOrientA.Ev1Mv2(molAleafPos1, molAleafPos0);
		tempOrientA.normalize();
		int molAOrienti =getOrientationIndex(tempOrientA);
		molAVec.ME(latticeOffset[molAOrienti]);
		
		// Molecule B
		molBVec.E(positionDefinition.position(moleculeB));

		Vector molBleafPos0 = moleculeB.getChildList().getAtom(0).getPosition();
	    Vector molBleafPos1 = moleculeB.getChildList().getAtom(1).getPosition();
	    	 
		tempOrientB. Ev1Mv2(molBleafPos1, molBleafPos0);
		tempOrientB.normalize();		
		int molBOrienti =getOrientationIndex(tempOrientB);
		molBVec.ME(latticeOffset[molBOrienti]);
				
		if(isReverseOrder){
			tempVec.Ev1Mv2(molAVec, molBVec);
		} else {
			tempVec.Ev1Mv2(molBVec, molAVec);	
		}
	
		coordinateDef.getBox().getBoundary().nearestImage(tempVec);
		int[] siteIndex = getSiteDisplacementIndex(tempVec);
		
		//System.out.println("tempVec: " + tempVec.toString());
		index[0] = siteIndex[0];
		index[1] = siteIndex[1];
		index[2] = siteIndex[2];
		
		if(isReverseOrder){
			index[3] = getOrientationIndex(tempOrientB);
			index[4] = getOrientationIndex(tempOrientA);
			
		} else{
			index[3] = getOrientationIndex(tempOrientA);
			index[4] = getOrientationIndex(tempOrientB);
				
		}
		//System.out.println(index[0]+"  "+index[1]+"  "+index[2]);
		return index;
	}
		
	public int[] getSiteDisplacementIndex(Vector siteDisplacement){
		int[] index = new int[3];
		index[0] = (int)Math.round(tempVec.getX(0)/lengthX) +  nCell;
		index[1] = (int)Math.round(tempVec.getX(1)/lengthY) + (nCell*2);
		index[2] = (int)Math.round(tempVec.getX(2)/lengthZ) +  nCell;
		
		return index; 
	}
	
	public void updateNewMoleculePair(int[] index){
		isNewPair[index[0]][index[1]][index[2]][index[3]][index[4]]=false;
		
	}
	
	public void resetNewMoleculePair(){
		for(int a=0; a<isNewPair.length; a++){
			for(int b=0; b<isNewPair[0].length; b++){
				for(int c=0; c<isNewPair[0][0].length; c++){
					for(int d=0; d<isNewPair[0][0][0].length; d++){
						for(int e=0; e<isNewPair[0][0][0][0].length; e++){
							isNewPair[a][b][c][d][e] = true;
						}	
					}	
				}	
			}	
		}
	}
	
	public int getOrientationIndex(Vector orientation){
		double x = orientation.getX(0);
		double y = orientation.getX(1);
		double z = orientation.getX(2);
		
		if(x > 0.0 && y > 0.0 && z < 0.0){
			return 0;
			
		} else if (y > 0.0 && z > 0.0){
			return 1;
			
		} else if (y < 0.0 && z > 0.0) {
			return 2;
			
		} else {
			return 3;
			
		}
	}
	
	public boolean getIsNewPair(int[] index) {
		return isNewPair[index[0]][index[1]][index[2]][index[3]][index[4]];
	}
	
	protected int[] index;
	protected CoordinateDefinitionNitrogen coordinateDef;
	protected IMoleculePositionDefinition positionDefinition;
	protected Vector tempVec, tempOrientA, tempOrientB, molAVec, molBVec;
	protected Vector[] orientation, latticeOffset;
	protected double lengthX, lengthY, lengthZ;
	protected double[][][] siteDisplacement;
	protected boolean[][][][][] isNewPair;
	protected int nCell;
}
