package etomica.models.nitrogen;

import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IVectorMutable;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.IAtomPositionDefinition;
import etomica.space.ISpace;


/**
 * 
 * 
 * @author taitan
 *
 */
public class FindPairMoleculeIndex {
	public FindPairMoleculeIndex(ISpace space, CoordinateDefinitionNitrogen coordinateDefination){
		this.coordinateDef = coordinateDefination;
		molecules = coordinateDefination.getBox().getMoleculeList();
		positionDefinition = new AtomPositionGeometricCenter(space);
		tempVec = space.makeVector();
		tempOrientA = space.makeVector();
		tempOrientB = space.makeVector();
		
		molAVec = space.makeVector();
		molBVec = space.makeVector();
		
		nCell = (int) Math.round(Math.pow((molecules.getMoleculeCount()/4), 1.0/3.0));
		halfUnitCellLength = coordinateDefination.getBox().getBoundary().getBoxSize().getX(0)/(2.0*nCell);
		
		int nSites = 2*nCell+1;
		siteDisplacement = new double[nSites][nSites][nSites];
		orientation = new IVectorMutable[4];
		
		for (int i=0; i<orientation.length; i++){
			orientation[i] = space.makeVector();
		}

		indices = new int [nSites][nSites][nSites][4][4][molecules.getMoleculeCount()][2];
		counter = new int [nSites][nSites][nSites][4][4];
		
		for (int a=0; a<counter.length; a++){
			for (int b=0; b<counter[0].length; b++){
				for (int c=0; c<counter[0][0].length; c++){
					for (int d=0; d<counter[0][0][0].length; d++){
						for (int e=0; e<counter[0][0][0][0].length; e++){
							counter[a][b][c][d][e] = 0;
						}	
					}	
				}	
			}	
		}
		
	}
	
	public void screenMolecules(){
		int numMol = molecules.getMoleculeCount();
		
		for (int i=0; i<numMol; i++){
			IMolecule moleculeA = molecules.getMolecule(i);
			molAVec.E(positionDefinition.position(moleculeA));
			
			IVectorMutable molAleafPos0 = moleculeA.getChildList().getAtom(0).getPosition();
    	   	IVectorMutable molAleafPos1 = moleculeA.getChildList().getAtom(1).getPosition();
    	 
			tempOrientA.Ev1Mv2(molAleafPos1, molAleafPos0);
			tempOrientA.normalize();
    	   	
			for(int j=i; j<numMol; j++){
				IMolecule moleculeB = molecules.getMolecule(j);
				molBVec.E(positionDefinition.position(moleculeB));
				tempVec.Ev1Mv2(molBVec, molAVec);
				coordinateDef.getBox().getBoundary().nearestImage(tempVec);
				
				IVectorMutable molBleafPos0 = moleculeB.getChildList().getAtom(0).getPosition();
	    	   	IVectorMutable molBleafPos1 = moleculeB.getChildList().getAtom(1).getPosition();
	    	 
				tempOrientB. Ev1Mv2(molBleafPos1, molBleafPos0);
				tempOrientB.normalize();
				
				int[] siteIndex = getSiteDisplacementIndex(tempVec);
				int x = siteIndex[0];
				int y = siteIndex[1];
				int z = siteIndex[2];
				int orientA = getOrientationIndex(tempOrientA);
				int orientB = getOrientationIndex(tempOrientB);
				
				//System.out.println(x+" "+y+ " " + z+ " "+orientA + " " + orientB);
				counter[x][y][z][orientA][orientB]++;
				indices[x][y][z][orientA][orientB][counter[x][y][z][orientA][orientB]][0] = i;
				indices[x][y][z][orientA][orientB][counter[x][y][z][orientA][orientB]][1] = j;
				//System.out.println(i +" " +j +" vector: " + tempVec.getX(0)+ " " + getOrientationIndex(tempOrient));
//				System.out.println(i +" " +j +" vector: " + siteIndex[0] +" " + siteIndex[1] +" " +siteIndex[2] 
//				              + " " + orientA+ " " + orientB);
			}
			//System.out.println("**************\n");
		}
		
	}
	
	public int[] getSiteDisplacementIndex(IVectorMutable siteDisplacement){
		int[] index = new int[3];
		index[0] = (int)Math.round(tempVec.getX(0)/halfUnitCellLength) + nCell;
		index[1] = (int)Math.round(tempVec.getX(1)/halfUnitCellLength) + nCell;
		index[2] = (int)Math.round(tempVec.getX(2)/halfUnitCellLength) + nCell;
		
		return index;
	}
	
	public int getOrientationIndex(IVectorMutable orientation){
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
	
	public int[][][][][][][] getIndices() {
		return indices;
	}

	public void setIndices(int[][][][][][][] indices) {
		this.indices = indices;
	}

	public int[][][][][] getCounter() {
		return counter;
	}

	public void setCounter(int[][][][][] counter) {
		this.counter = counter;
	}

	protected int[][][][][][][] indices;
	protected int[][][][][] counter;
	protected CoordinateDefinitionNitrogen coordinateDef;
	protected IMoleculeList molecules;
	protected IAtomPositionDefinition positionDefinition;
	protected IVectorMutable tempVec, tempOrientA, tempOrientB, molAVec, molBVec;
	protected IVectorMutable[] orientation;
	protected double halfUnitCellLength;
	protected double[][][] siteDisplacement;
	protected int nCell;
}
