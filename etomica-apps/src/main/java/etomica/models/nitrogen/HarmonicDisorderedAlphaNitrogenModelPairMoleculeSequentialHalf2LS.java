/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeChildAtomAction;
import etomica.atom.IMolecule;
import etomica.api.ISpecies;
import etomica.space.Vector;
import etomica.atom.MoleculePair;
import etomica.box.Box;
import etomica.data.types.DataTensor;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.normalmode.BasisBigCell;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;

/**
 * Disordered alpha-phase
 * 
 * This class is created to take care of Java out of memory problem when creating
 *  a 2D array that is gigantic [~10000][~10000]
 *  
 * The class is write a 3 x (nA * dof) Matrix to file while looping through the molecules
 * This class also applied the molecules pair identification algorithm to speed up the
 *  matrix construction. 
 *  
 * Additional: For larger system size, the matrix can be broken down into smaller matrices 
 *              and then can be summed up later, and it is called the Half2
 * 
 * with lattice sum
 * 
 * @author Tai Boon Tan
 *
 */
public class HarmonicDisorderedAlphaNitrogenModelPairMoleculeSequentialHalf2LS extends Simulation{

	
	public HarmonicDisorderedAlphaNitrogenModelPairMoleculeSequentialHalf2LS(Space space, int numMolecule, double density, double rC) {
		super(space);
		this.space = space;
		
		int nCell = (int) Math.round(Math.pow((numMolecule/4), 1.0/3.0));
		double unitCellLength = Math.pow(numMolecule/density, 1.0/3.0)/nCell;//5.661;
//		System.out.println("a: " + unitCellLength);
//		System.out.println("nCell: " + nCell);
		
		potentialMaster = new PotentialMaster();
				
		int division = 2;
		Basis basisFCC = new BasisCubicFcc();
		Basis basis = new BasisBigCell(space, basisFCC, new int[]{nCell/division, nCell/division, nCell/division});
		
		ConformationNitrogen conformation = new ConformationNitrogen(space);
		species = new SpeciesN2(space);
		species.setConformation(conformation);
		addSpecies(species);
		
		box = new Box(space);
		addBox(box);
		box.setNMolecules(species, numMolecule);		
		
		int [] nCells = new int[]{division,division,division};
		Boundary boundary = new BoundaryRectangularPeriodic(space);
		boundary.setBoxSize(space.makeVector(new double[]{nCell*unitCellLength,nCell*unitCellLength,nCell*unitCellLength}));
		Primitive primitive = new PrimitiveCubic(space, (nCell/division)*unitCellLength);
	
		coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordinateDef.setIsAlpha();
		coordinateDef.setOrientationVectorAlpha(space);
		coordinateDef.initializeCoordinates(nCells);
		
		box.setBoundary(boundary);
//		System.out.println("Truncation Radius (" + rCScale +" Box Length): " + rC);
		
		potential = new P2Nitrogen(space, rC);
		potential.setEnablePBC(false);
		potential.setBox(box);

		potentialMaster.addPotential(potential, new ISpecies[]{species, species});
		
		int nSites = 2*nCell+1;
		pairMatrix = new double[nSites][nSites][nSites][4][4][3][3];
		isFoundReverse = new boolean[nSites][nSites][nSites][4][4]; //default to false

		findPair = new FindPairMoleculeIndex(space, coordinateDef);
		
		translateBy = new AtomActionTranslateBy(coordinateDef.getPrimitive().getSpace());
        atomGroupActionTranslate = new MoleculeChildAtomAction(translateBy); 
		lsPosition = coordinateDef.getPrimitive().getSpace().makeVector();
        
		xVecBox = coordinateDef.getBox().getBoundary().getBoxSize().getX(0);
		yVecBox = coordinateDef.getBox().getBoundary().getBoxSize().getX(1);
		zVecBox = coordinateDef.getBox().getBoundary().getBoxSize().getX(2); 
		
		double rX = coordinateDef.getBox().getBoundary().getBoxSize().getX(0);
		this.nLayer = (int)Math.round(rC/rX + 0.5);
		if (this.nLayer < 1){
			throw new RuntimeException("<HarmonicDisorderedAlphaNitrogenModelPairMoleculeSequentialHalf2LS> nLayer is "+ this.nLayer);
		}
		
		workArray = new double[3][3*numMolecule];
//		System.out.println("rX: " + rX);
//		System.out.println("nLayer: " + nLayer);
	}
	
	public double[][] get2ndDerivative(int molec0){
	
		int numMolecule = box.getMoleculeList().getMoleculeCount();
		int dofTrans = 3;	
//		workArray = new double[3][3*numMolecule];
		
		DataTensor transTensor = new DataTensor(space);
		MoleculePair pair = new MoleculePair();	
		
		for (int i=0; i<3; i++){
			for (int j=0; j<3*numMolecule; j++){
				if((j/3) == molec0){
					workArray[i][j] = 0.0;
				}else {
					workArray[i][j] = Double.NaN;
				}
				
			}
		}
		/*
		 *	Constructing the upper diagonal of the matrix
		 *	(Skipping the molec1 == molec2) 
		 */
		IMolecule molecule0 = coordinateDef.getBox().getMoleculeList().getMolecule(molec0);
		pair.atom0 = molecule0;
		
		boolean isReverseOrder = false;
		for(int molec1=0; molec1<numMolecule; molec1++){
			if(molec0 == molec1) continue;
			
			if(molec0 > molec1){
				isReverseOrder = true;
			} else {
				isReverseOrder = false;
			}
			
			
			// Analytical calculation for 3x3 Translational second Derivative
			IMolecule molecule1 = coordinateDef.getBox().getMoleculeList().getMolecule(molec1);
			pair.atom1 = molecule1;
			
			int[] index = findPair.getPairMoleculesIndex(pair.atom0, pair.atom1, isReverseOrder);
			boolean isNewPair = findPair.getIsNewPair(index);
			
			if(isReverseOrder && isNewPair){
				isFoundReverse[index[0]][index[1]][index[2]][index[3]][index[4]] = true;
			}
			
			if(isNewPair){
				
				transTensor.E(0.0);
				
				//do Lattice sum
				for(int x=-nLayer; x<=nLayer; x++){
					for(int y=-nLayer; y<=nLayer; y++){
						for(int z=-nLayer; z<=nLayer; z++){
							lsPosition.E(new double[]{x*xVecBox, y*yVecBox, z*zVecBox});
							translateBy.setTranslationVector(lsPosition);
							atomGroupActionTranslate.actionPerformed(molecule1);
		
							transTensor.PE(potential.secondDerivative(pair));
														
							lsPosition.TE(-1);
							translateBy.setTranslationVector(lsPosition);
							atomGroupActionTranslate.actionPerformed(molecule1);
						}	
					}	
				}
				
				for(int i=0; i<3; i++){
					for(int j=0; j<3; j++){
						workArray[i][molec1*dofTrans + j] = transTensor.x.component(i, j);
						pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j] = workArray[i][molec1*dofTrans + j];
					}
				}
				
				findPair.updateNewMoleculePair(index);
			} else {
				
				if(isFoundReverse[index[0]][index[1]][index[2]][index[3]][index[4]]==true){
					isReverseOrder = !isReverseOrder; 
				}
				
				if(isReverseOrder){
					for(int i=0; i<3; i++){
						for(int j=0; j<3; j++){
							workArray[i][molec1*dofTrans + j] 
							                            = pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][j][i];
						}
					}
				} else {
					for(int i=0; i<3; i++){
						for(int j=0; j<3; j++){
							workArray[i][molec1*dofTrans + j] 
							                            = pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j];
						}
					}
					
				}
			}
		}
			
      	/*
      	 *  SELF-TERM
      	 *  
      	 *  The diagonal-element block are all the same
      	 *  The 3x3 translation block is found by summing all the interactions with molecules in the system
      	 *   however, the cross and rotation block have be determined numerically 
      	 */

		int[] index = findPair.getPairMoleculesIndex(molecule0, molecule0, false);
		boolean isNewPair = findPair.getIsNewPair(index);
		
		if(isNewPair){
			for(int molec1=0; molec1<numMolecule; molec1++){
    			if(molec0==molec1) continue; // we might double sum the elements in array[a][a] if we don't skip the pair
    			for(int i=0; i<3; i++){
    				for (int j=0; j<3; j++){
    					workArray[i][molec0*dofTrans + j] -= workArray[i][molec1*dofTrans + j] ;         				
	    				
    				}
    			}
    		}
			
			for(int i=0; i<3; i++){
   				for (int j=0; j<3; j++){
   					pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j] = workArray[i][molec0*dofTrans + j];
   				}
			}
			
    	
		} else {
				
			for(int i=0; i<3; i++){
   				for (int j=0; j<3; j++){
   					workArray[i][molec0*dofTrans + j]= pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j];
   				}
			}
			
		}
		
		for (int i=0; i<3; i++){
			for (int j=0; j<3*numMolecule; j++){
				if (Double.isNaN(workArray[i][j])){
					throw new RuntimeException("BUSTED !!! element["+i+"]["+j+"] is Nan");
				}
			}
		}
		return workArray;
		
	}
	
	public void constructHessianMatrix(int nCell){

		int numMolecules = nCell*nCell*nCell*4;

		for (int iMol=0; iMol<numMolecules; iMol++){
			double[][]array = get2ndDerivative(iMol);
			
			for (int iRow=0; iRow< array.length; iRow++){
				for (int jCol=0; jCol< array[0].length; jCol++){
					double value = array[iRow][jCol];
//						if(Math.abs(value) < 1e-6){
//							value = 0.0;
//						}
					System.out.print(value+ " ");
				}
				System.out.println("");
			}
		}
	}
	
	public static void main (String[] args){
		
		int nC=3;
		double density = 0.022;
		double rC = 100;
		if(args.length > 0){
			nC = Integer.parseInt(args[0]);
		}
		if(args.length > 1){
			density = Double.parseDouble(args[1]);
		}
		if(args.length > 2){
			rC = Double.parseDouble(args[2]);
		}
		
		int numMolecule =nC*nC*nC*4;

		HarmonicDisorderedAlphaNitrogenModelPairMoleculeSequentialHalf2LS test = new HarmonicDisorderedAlphaNitrogenModelPairMoleculeSequentialHalf2LS(Space3D.getInstance(3), numMolecule, density, rC);
		test.constructHessianMatrix(nC);
		
	}
	
	
	protected Box box;
	protected Space space;
	protected P2Nitrogen potential;
	protected CoordinateDefinitionNitrogen coordinateDef;
	protected PotentialMaster potentialMaster;
	protected double[][][][][][][] pairMatrix;
	protected FindPairMoleculeIndex findPair;
	protected SpeciesN2 species;
	protected AtomActionTranslateBy translateBy;
	protected MoleculeChildAtomAction atomGroupActionTranslate;
	protected Vector lsPosition;
	protected double xVecBox, yVecBox, zVecBox, rC;
	protected int nLayer;
	protected boolean[][][][][] isFoundReverse;
	protected double[][] workArray;
	
	private static final long serialVersionUID = 1L;
}
