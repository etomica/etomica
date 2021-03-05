/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.data.types.DataTensor;
import etomica.graphics.DisplayBox;
import etomica.graphics.SimulationGraphic;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculePair;
import etomica.normalmode.BasisBigCell;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesGeneral;
import etomica.units.Pixel;

/**
 * This class is created to take care of Java out of memory problem when creating
 *  a 2D array that is gigantic [~10000][~10000]
 *  
 * The class is write a 5 x (nA * dof) Matrix to file while looping through the molecules
 * This class also applied the molecules pair identification algorithm to speed up the
 *  matrix construction. 
 *  
 * Additional: For larger system size, the matrix can be broken down into smaller matrices 
 *              and then can be summed up later, and it is called the Half2
 * 
 * @author Tai Boon Tan
 *
 */
public class HarmonicAlphaNitrogenModelPairMoleculeSequentialHalf2 extends Simulation{

	
	public HarmonicAlphaNitrogenModelPairMoleculeSequentialHalf2(Space space, int numMolecule, double density) {
		super(space);
		this.space = space;

        species = SpeciesN2.create(false);
        addSpecies(species);

		int nCell = (int) Math.round(Math.pow((numMolecule / 4), 1.0 / 3.0));
		double unitCellLength = Math.pow(numMolecule / density, 1.0 / 3.0) / nCell;//5.661;
//		System.out.println("a: " + unitCellLength);
//		System.out.println("nCell: " + nCell);

		potentialMaster = new PotentialMaster();

		int division = 2;
        Basis basisFCC = new BasisCubicFcc();
        Basis basis = new BasisBigCell(space, basisFCC, new int[]{nCell / division, nCell / division, nCell / division});

        Boundary boundary = new BoundaryRectangularPeriodic(space, nCell * unitCellLength);
        box = this.makeBox(boundary);
        box.setNMolecules(species, numMolecule);

        int[] nCells = new int[]{division, division, division};
        Primitive primitive = new PrimitiveCubic(space, (nCell / division) * unitCellLength);

        coordinateDef = new CoordinateDefinitionNitrogen(getSpeciesManager(), box, primitive, basis, space);
        coordinateDef.setIsAlpha();
        coordinateDef.setOrientationVectorAlpha(space);
        coordinateDef.initializeCoordinates(nCells);
        double rCScale = 0.475;
        double rC = box.getBoundary().getBoxSize().getX(0) * rCScale;
//		System.out.println("Truncation Radius (" + rCScale +" Box Length): " + rC);

        potential = new P2Nitrogen(space, rC);
        potential.setBox(box);

		potentialMaster.addPotential(potential, new ISpecies[]{species, species});

		int nSites = 2 * nCell + 1;
		pairMatrix = new double[nSites][nSites][nSites][4][4][5][5];

//		cm2ndD = new CalcNumerical2ndDerivativeNitrogen(box, potential, coordinateDef);

		cAN2nD = new CalcHalfAnalyticHalfNumeric2ndDerivativeNitrogen(space, box, potential, coordinateDef, true);
		cA2nD = new CalcAnalytical2ndDerivativeNitrogen(space, box, potential, coordinateDef);
		findPair = new FindPairMoleculeIndex(space, coordinateDef);
	}
	
	public double[][] get2ndDerivative(int molec0){
	
		DataTensor transTensor = new DataTensor(space);
		MoleculePair pair = new MoleculePair();
	
		int numBasis = coordinateDef.getBasis().getScaledCoordinates().length;
		int numBasisCell = coordinateDef.getBasisCells().length;
		int numMolecules = numBasis * numBasisCell;
		int dofPerMol = coordinateDef.getCoordinateDim()/numBasis;
		double[][] array = new double[dofPerMol][numBasisCell*coordinateDef.getCoordinateDim()];
			
		/*
		 *	Constructing the upper diagonal of the matrix
		 *	(Skipping the molec1 == molec2) 
		 */
		IMolecule molecule0 = coordinateDef.getBox().getMoleculeList().get(molec0);
		pair.mol0 = molecule0;
			
		boolean isReverseOrder = false;
		for(int molec1=0; molec1<numMolecules; molec1++){
			if(molec0 == molec1) continue;
			/*
			 * working within the 5x5 Matrix
			 */
			// Analytical calculation for 3x3 Translational second Derivative
			
			if(molec0 > molec1){
				isReverseOrder = true;
			} else {
				isReverseOrder = false;
			}
			
			pair.mol1 = coordinateDef.getBox().getMoleculeList().get(molec1);
		
			int[] index = findPair.getPairMoleculesIndex(pair.mol0, pair.mol1, isReverseOrder);
			boolean isNewPair = findPair.getIsNewPair(index);
				
			if(isNewPair){
				double[][] a = cA2nD.d2phi_du2(new int[]{molec1,molec0});
				for(int i=0; i<dofPerMol; i++){
					for(int j=0; j<dofPerMol; j++){
						array[i][molec1*dofPerMol + j] = a[i][j];
						pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j] = array[i][molec1*dofPerMol + j];
					}
				}
									
				findPair.updateNewMoleculePair(index);
					
			} else {
				
				if(isReverseOrder){
					for(int i=0; i<3; i++){
						for(int j=0; j<3; j++){
							array[i][molec1*dofPerMol + j] 
							                            = pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][j][i];
						}
					}
					// Numerical calculation for the Cross (trans and rotation) and rotation second Derivative
					for(int i=0; i<dofPerMol; i++){
						for(int j=0; j<dofPerMol; j++){
							if(i<3 && j<3) continue;
							array[i][molec1*dofPerMol + j] = 
								pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][j][i];
						}
					}
				} else {
					for(int i=0; i<3; i++){
						for(int j=0; j<3; j++){
							array[i][molec1*dofPerMol + j] 
							                            = pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j];
						}
					}
					// Numerical calculation for the Cross (trans and rotation) and rotation second Derivative
					for(int i=0; i<dofPerMol; i++){
						for(int j=0; j<dofPerMol; j++){
							if(i<3 && j<3) continue;
							array[i][molec1*dofPerMol + j] = 
								pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j];
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
			for(int molec1=0; molec1<numMolecules; molec1++){
    			if(molec0==molec1) continue; // we might double sum the elements in array[a][a] if we don't skip the pair
	    			
    			for(int i=0; i<3; i++){
    				for (int j=0; j<3; j++){
    					array[i][molec0*dofPerMol + j] -= array[i][molec1*dofPerMol + j] ;         				
	    				
    				}
    			}
    		}
				
			for(int i=0; i<3; i++){
   				for (int j=0; j<3; j++){
   					pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j] = array[i][molec0*dofPerMol + j];
   				}
			}
    			
		} else {
				
			for(int i=0; i<3; i++){
   				for (int j=0; j<3; j++){
   					array[i][molec0*dofPerMol + j]= pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j];
   				}
			}
			
		}
			
		// Numerical calculation for the Cross (trans and rotation) and rotation second Derivative
			
		if(isNewPair){
			for(int i=0; i<dofPerMol; i++){
				for(int j=0; j<dofPerMol; j++){
					if(i<3 && j<3) continue;
					// j i because it got switched molecule A and molecule B
					array[i][molec0*dofPerMol + j] = cAN2nD.d2phi_du2(new int[]{molec0,molec0}, new int[]{j,i});
					pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j] = array[i][molec0*dofPerMol + j];
				}    		
	    	}
			findPair.updateNewMoleculePair(index);
				
		} else {
			for(int i=0; i<dofPerMol; i++){
				for(int j=0; j<dofPerMol; j++){
					if(i<3 && j<3) continue;
					array[i][molec0*dofPerMol + j] =  pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j];
				}    		
	   		}
			
		}
		
		return array;
		
	}
	
	public void constructHessianMatrix(int nCell){
		
		int interval = nCell*2;
		int numBasis = coordinateDef.getBasis().getScaledCoordinates().length;
		int numBasisCell = coordinateDef.getBasisCells().length;
		
		int numMolecules = numBasis * numBasisCell;
		
		int dofPerMol = coordinateDef.getCoordinateDim()/numBasis;	
		double[][][] array = new double[interval][dofPerMol][numBasisCell*coordinateDef.getCoordinateDim()];
			
		for (int iMol=0; iMol<numMolecules; iMol+=interval){
			for(int i=0; i<interval; i++){
				array[i] = get2ndDerivative(iMol+i);
			}
			for(int i=0; i<interval; i++){
				for (int iRow=0; iRow< array[0].length; iRow++){
					for (int jCol=0; jCol< array[0][0].length; jCol++){
						double value = array[i][iRow][jCol];
//							if(Math.abs(value) < 1e-6){
//								value = 0.0;
//							}
						System.out.print(value+ " ");
					}
					System.out.println("");
				}
			}
		}
			
	
	}
	
	public static void main (String[] args){
		
		int nC=4;
		double density = 0.0240;
		if(args.length > 0){
			nC = Integer.parseInt(args[0]);
		}
		if(args.length > 1){
			density = Double.parseDouble(args[1]);
		}
		
		int numMolecule =nC*nC*nC*4;

		HarmonicAlphaNitrogenModelPairMoleculeSequentialHalf2 test = new HarmonicAlphaNitrogenModelPairMoleculeSequentialHalf2(Space3D.getInstance(3), numMolecule, density);
		test.constructHessianMatrix(nC);
	
		if(false){
			SimulationGraphic simGraphic = new SimulationGraphic(test, SimulationGraphic.TABBED_PANE);
			simGraphic.add(new DisplayBox(test, test.box));
			simGraphic.getDisplayBox(test.box).setPixelUnit(new Pixel(10));
			
			DiameterHashByType diameter = new DiameterHashByType();
			diameter.setDiameter(test.species.getTypeByName("N"), 3.1);
			diameter.setDiameter(test.species.getTypeByName("P"), 0.0);

			simGraphic.getDisplayBox(test.box).setDiameterHash(diameter);
			
	
		    simGraphic.makeAndDisplayFrame("Alpha Crystal Structure");
		    return;
		}
		
	}
	
	
	protected Box box;
	protected Space space;
	protected P2Nitrogen potential;
	protected CoordinateDefinitionNitrogen coordinateDef;
	protected PotentialMaster potentialMaster;
	protected double[][][][][][][] pairMatrix;
//	protected CalcNumerical2ndDerivativeNitrogen cm2ndD;
	protected CalcHalfAnalyticHalfNumeric2ndDerivativeNitrogen cAN2nD;
	protected CalcAnalytical2ndDerivativeNitrogen cA2nD;
	protected FindPairMoleculeIndex findPair;
	protected SpeciesGeneral species;
	private static final long serialVersionUID = 1L;
}
