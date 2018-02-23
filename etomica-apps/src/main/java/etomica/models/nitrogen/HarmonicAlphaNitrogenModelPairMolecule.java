/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.box.Box;
import etomica.data.types.DataTensor;
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

import java.io.FileWriter;
import java.io.IOException;



/**
 * 
 * @author Tai Boon Tan
 *
 */
public class HarmonicAlphaNitrogenModelPairMolecule extends Simulation{

	
	public HarmonicAlphaNitrogenModelPairMolecule(Space space, int numMolecule, double density) {
        super(space);
        this.space = space;


        int nCell = (int) Math.round(Math.pow((numMolecule / 4), 1.0 / 3.0));
        double unitCellLength = Math.pow(numMolecule / density, 1.0 / 3.0) / nCell;//5.661;
        System.out.println("a: " + unitCellLength);
        System.out.println("nCell: " + nCell);

        potentialMaster = new PotentialMaster();

        Basis basisFCC = new BasisCubicFcc();
        Basis basis = new BasisBigCell(space, basisFCC, new int[]{nCell, nCell, nCell});

        ConformationNitrogen conformation = new ConformationNitrogen(space);
        SpeciesN2 species = new SpeciesN2(space);
        species.setConformation(conformation);
        addSpecies(species);

        box = this.makeBox();
        box.setNMolecules(species, numMolecule);

        int[] nCells = new int[]{1, 1, 1};
        Boundary boundary = new BoundaryRectangularPeriodic(space, nCell * unitCellLength);
        Primitive primitive = new PrimitiveCubic(space, nCell * unitCellLength);

        coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
        coordinateDef.setIsAlpha();
        coordinateDef.setOrientationVectorAlpha(space);
        coordinateDef.initializeCoordinates(nCells);

        box.setBoundary(boundary);
        double rCScale = 0.475;
        double rC = box.getBoundary().getBoxSize().getX(0) * rCScale;
        System.out.println("Truncation Radius (" + rCScale + " Box Length): " + rC);

        potential = new P2Nitrogen(space, rC);
        potential.setBox(box);

        potentialMaster.addPotential(potential, new ISpecies[]{species, species});

        int nSites = 2 * nCell + 1;
        pairMatrix = new double[nSites][nSites][nSites][4][4][5][5];

    }
	
	public double[][] get2ndDerivative(){
		
		double[][] array = new double[coordinateDef.getCoordinateDim()][coordinateDef.getCoordinateDim()];
		DataTensor transTensor = new DataTensor(space);
		MoleculePair pair = new MoleculePair();
	
		int numMolecule = box.getMoleculeList().size();
		int dofPerMol = coordinateDef.getCoordinateDim()/numMolecule;
		
		CalcNumerical2ndDerivativeNitrogen cm2ndD = new CalcNumerical2ndDerivativeNitrogen(box, potential, coordinateDef);
		FindPairMoleculeIndex findPair = new FindPairMoleculeIndex(space, coordinateDef);
		
		/*
		 *	Constructing the upper diagonal of the matrix
		 *	(Skipping the molec1 == molec2) 
		 */
		for(int molec0=0; molec0<numMolecule; molec0++){
			pair.mol0 = box.getMoleculeList().get(molec0);
			
			for(int molec1=molec0; molec1<numMolecule; molec1++){
				if(molec0 == molec1) continue;
				/*
				 * working within the 5x5 Matrix
				 */
				// Analytical calculation for 3x3 Translational second Derivative
				pair.mol1 = box.getMoleculeList().get(molec1);
		
				int[] index = findPair.getPairMoleculesIndex(pair.mol0, pair.mol1, false);
				boolean isNewPair = findPair.getIsNewPair(index);
				
				if(isNewPair){
					transTensor.E(potential.secondDerivative(pair));
					for(int i=0; i<3; i++){
						for(int j=0; j<3; j++){
							array[molec0*dofPerMol + i][molec1*dofPerMol + j] = transTensor.x.component(i, j);
							pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j] = array[molec0*dofPerMol + i][molec1*dofPerMol + j];
						}
					}
					// Numerical calculation for the Cross (trans and rotation) and rotation second Derivative
					for(int i=0; i<dofPerMol; i++){
						for(int j=0; j<dofPerMol; j++){
							if(i<3 && j<3) continue;
							// j i because it got switched molecule A and molecule B
							array[molec0*dofPerMol + i][molec1*dofPerMol + j] = cm2ndD.d2phi_du2(new int[]{molec0,molec1}, new int[]{j,i});
							pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j] = array[molec0*dofPerMol + i][molec1*dofPerMol + j];
						}
					}
					
					findPair.updateNewMoleculePair(index);
					
				} else {
					for(int i=0; i<3; i++){
						for(int j=0; j<3; j++){
							array[molec0*dofPerMol + i][molec1*dofPerMol + j] 
							                            = pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j];
						}
					}
					// Numerical calculation for the Cross (trans and rotation) and rotation second Derivative
					for(int i=0; i<dofPerMol; i++){
						for(int j=0; j<dofPerMol; j++){
							if(i<3 && j<3) continue;
							array[molec0*dofPerMol + i][molec1*dofPerMol + j] = 
								pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j];
						}
					}
				}
				
			}
		}
		
		/*
		 *  Assigning the lower diagonal of the matrix
		 *  (Skipping the molec1 == molec2) 
		 */
		for(int molec1=0; molec1<numMolecule; molec1++){
			for(int molec0=molec1; molec0<numMolecule; molec0++){
				
				if(molec0 == molec1) continue;
				for(int i=0; i<dofPerMol; i++){
					for (int j=0; j<dofPerMol; j++){
						array[molec0*dofPerMol + i][molec1*dofPerMol + j] = array[molec1*dofPerMol + j][molec0*dofPerMol + i];
						
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

		for(int molec0=0; molec0<numMolecule; molec0++){
			IMolecule molecule0 = coordinateDef.getBox().getMoleculeList().get(molec0);
			int[] index = findPair.getPairMoleculesIndex(molecule0, molecule0, false);
			boolean isNewPair = findPair.getIsNewPair(index);
			
			if(isNewPair){
				for(int molec1=0; molec1<numMolecule; molec1++){
	    			if(molec0==molec1) continue; // we might double sum the elements in array[a][a] if we don't skip the pair
	    			
	    			for(int i=0; i<3; i++){
	    				for (int j=0; j<3; j++){
	    					array[molec0*dofPerMol + i][molec0*dofPerMol + j] -= array[molec0*dofPerMol + i][molec1*dofPerMol + j] ;         				
	    				
	    				}
	    			}
	    		}
				
				for(int i=0; i<3; i++){
    				for (int j=0; j<3; j++){
    					pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j] = array[molec0*dofPerMol + i][molec0*dofPerMol + j];
    				}
				}
    			
			} else {
				
				for(int i=0; i<3; i++){
    				for (int j=0; j<3; j++){
    					array[molec0*dofPerMol + i][molec0*dofPerMol + j]= pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j];
    				}
				}
			
			}
			
    	}
				
		for(int molec0=0; molec0<numMolecule; molec0++){
			IMolecule molecule0 = coordinateDef.getBox().getMoleculeList().get(molec0);
			int[] index = findPair.getPairMoleculesIndex(molecule0, molecule0, false);
			boolean isNewPair = findPair.getIsNewPair(index);
			
			// Numerical calculation for the Cross (trans and rotation) and rotation second Derivative
			
			if(isNewPair){
				for(int i=0; i<dofPerMol; i++){
					for(int j=0; j<dofPerMol; j++){
						if(i<3 && j<3) continue;
						// j i because it got switched molecule A and molecule B
						array[molec0*dofPerMol + i][molec0*dofPerMol + j] = cm2ndD.d2phi_du2(new int[]{molec0,molec0}, new int[]{j,i});
						pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j] = array[molec0*dofPerMol + i][molec0*dofPerMol + j];
					}    		
	    		}
				findPair.updateNewMoleculePair(index);
				
			} else {
				for(int i=0; i<dofPerMol; i++){
					for(int j=0; j<dofPerMol; j++){
						if(i<3 && j<3) continue;
						array[molec0*dofPerMol + i][molec0*dofPerMol + j] =  pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j];
					}    		
	    		}
				
			}
    	}
		
		return array;
		
	}
	
	public static void main (String[] args){
		
		int nC=2;
		if(args.length > 0){
			nC = Integer.parseInt(args[0]);
		}
		int numMolecule =nC*nC*nC*4;
		double density = 0.025;
		HarmonicAlphaNitrogenModelPairMolecule test = new HarmonicAlphaNitrogenModelPairMolecule(Space3D.getInstance(3), numMolecule, density);
	
		long startTime = System.currentTimeMillis();
		
		double[] newU = new double[test.coordinateDef.getCoordinateDim()];
		double[][]testArray = test.get2ndDerivative();
		
		String fname = new String ("alpha"+numMolecule+"_2ndDer_d"+density+"_newPair");
		try {
			FileWriter fileWriter = new FileWriter(fname);
			
			for (int i=0; i<newU.length; i++){
				for (int j=0; j<newU.length; j++){
					double value = testArray[i][j];
//					if(Math.abs(value) < 1e-6){
//						value = 0.0;
//					}
					fileWriter.write(value+ " ");
				}
				fileWriter.write("\n");
			}
			fileWriter.close();
			
		} catch (IOException e) {
			
		}
		
		long endTime = System.currentTimeMillis();
		System.out.println("Time taken (s): " + (endTime-startTime)/1000);
	
	}
	
	
	protected Box box;
	protected Space space;
	protected P2Nitrogen potential;
	protected CoordinateDefinitionNitrogen coordinateDef;
	protected PotentialMaster potentialMaster;
	protected double[][][][][][][] pairMatrix;
	private static final long serialVersionUID = 1L;
}
