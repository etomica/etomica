/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.molecule.MoleculePositionGeometricCenter;
import etomica.space.Space;
import etomica.space.Vector;


/**
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
public class FindPairMoleculeIndex {
	public FindPairMoleculeIndex(Space space, CoordinateDefinitionNitrogen coordinateDefination){
		this.coordinateDef = coordinateDefination;
		positionDefinition = new MoleculePositionGeometricCenter(space);
		tempVec = space.makeVector();
		tempOrientA = space.makeVector();
		tempOrientB = space.makeVector();
		
		molAVec = space.makeVector();
		molBVec = space.makeVector();
		
		int numMolecule = coordinateDefination.getBox().getMoleculeList().size();
		nCell = (int) Math.round(Math.pow((numMolecule/4), 1.0/3.0));
		if(nCell > 20){
			throw new RuntimeException("<FindPairMoleculeIndex> nCell is greater than 20!!! " +
					"YOU ARE CRASHING JAVA MEMORY!! Live long and prosper!");
		}
		halfUnitCellLength = coordinateDefination.getBox().getBoundary().getBoxSize().getX(0)/(2.0*nCell);
		
		int nSites = 2*nCell+1;
		siteDisplacement = new double[nSites][nSites][nSites];

		index = new int[5];
		isNewPair = new boolean[nSites][nSites][nSites][4][4];
		
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
	
		molAVec.E(positionDefinition.position(moleculeA));
			
		Vector molAleafPos0 = moleculeA.getChildList().get(0).getPosition();
    	Vector molAleafPos1 = moleculeA.getChildList().get(1).getPosition();
    	 
		tempOrientA.Ev1Mv2(molAleafPos1, molAleafPos0);
		tempOrientA.normalize();
    	   	
		molBVec.E(positionDefinition.position(moleculeB));
		
		if(isReverseOrder){
			tempVec.Ev1Mv2(molAVec, molBVec);
		} else {
			tempVec.Ev1Mv2(molBVec, molAVec);	
		}
		coordinateDef.getBox().getBoundary().nearestImage(tempVec);
				
		Vector molBleafPos0 = moleculeB.getChildList().get(0).getPosition();
	    Vector molBleafPos1 = moleculeB.getChildList().get(1).getPosition();
	    	 
		tempOrientB. Ev1Mv2(molBleafPos1, molBleafPos0);
		tempOrientB.normalize();

		getSiteDisplacementIndex();
		
		if(isReverseOrder){
			index[3] = getOrientationIndex(tempOrientB);
			index[4] = getOrientationIndex(tempOrientA);
			
		} else{
			index[3] = getOrientationIndex(tempOrientA);
			index[4] = getOrientationIndex(tempOrientB);
				
		}
		return index;
	}
	
	public int[] getSiteDisplacementIndex(){
		index[0] = (int)Math.round(tempVec.getX(0)/halfUnitCellLength) + nCell;
		index[1] = (int)Math.round(tempVec.getX(1)/halfUnitCellLength) + nCell;
		index[2] = (int)Math.round(tempVec.getX(2)/halfUnitCellLength) + nCell;
		
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
		
		if(x > 0.0 && y > 0.0 && z > 0.0){
			return 0;
			
		} else if (y > 0.0 && z < 0.0){
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
	protected double halfUnitCellLength;
	protected double[][][] siteDisplacement;
	protected boolean[][][][][] isNewPair;
	protected int nCell;
}
