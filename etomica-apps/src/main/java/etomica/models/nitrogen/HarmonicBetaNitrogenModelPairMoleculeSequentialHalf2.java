/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.box.Box;
import etomica.data.types.DataTensor;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculePair;
import etomica.normalmode.BasisBigCell;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.units.Degree;

/**
 * 
 * Beta-phase Nitrogen Harmonic Approximation
 * with only 3 translational degrees of freedom
 * 
 * @author Tai Boon Tan
 *
 */
public class HarmonicBetaNitrogenModelPairMoleculeSequentialHalf2 extends Simulation{
	
	public HarmonicBetaNitrogenModelPairMoleculeSequentialHalf2(Space space, int numMolecule, String densityStr) {
        super(space);
        this.space = space;

        potentialMaster = new PotentialMaster();

        double density = Double.parseDouble(densityStr);
        int division = 2;
        double ratio = 1.631;
        double aDim = Math.pow(4.0 / (Math.sqrt(3.0) * ratio * density), 1.0 / 3.0);
        double cDim = aDim * ratio;
        //System.out.println("\naDim: " + aDim + " ;cDim: " + cDim);
        int nCell = (int) Math.pow(numMolecule / 1.999999999, 1.0 / 3.0);

        Basis basisHCP = new BasisHcp();
        BasisBigCell basis = new BasisBigCell(space, basisHCP, new int[]{nCell / division, nCell / division, nCell / division});

        Vector[] boxDim = new Vector[3];
        boxDim[0] = space.makeVector(new double[]{nCell * aDim, 0, 0});
        boxDim[1] = space.makeVector(new double[]{-nCell * aDim * Math.cos(Degree.UNIT.toSim(60)), nCell * aDim * Math.sin(Degree.UNIT.toSim(60)), 0});
        boxDim[2] = space.makeVector(new double[]{0, 0, nCell * cDim});

        int[] nCells = new int[]{division, division, division};
        Boundary boundary = new BoundaryDeformablePeriodic(space, boxDim);
        Primitive primitive = new PrimitiveHexagonal(space, (nCell / division) * aDim, (nCell / division) * cDim);

        SpeciesN2 species = new SpeciesN2(space);
        addSpecies(species);

        box = this.makeBox();
        box.setNMolecules(species, numMolecule);

        coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
        coordinateDef.setIsBeta();
        coordinateDef.setOrientationVectorBeta(space);
        coordinateDef.initializeCoordinates(nCells);
//		System.out.println("density: " + density);
//		System.exit(1);
        double[] u = new double[20];
        if (true) {
            BetaPhaseLatticeParameter parameters = new BetaPhaseLatticeParameter();
            double[][] param = parameters.getParameter(density);

            int kParam = 0;
            for (int i = 0; i < param.length; i++) {
                for (int j = 0; j < param[0].length; j++) {
                    u[kParam] = param[i][j];
                    kParam++;
                }
            }

            int numDOF = coordinateDef.getCoordinateDim();
            double[] newU = new double[numDOF];
            if (true) {
                for (int j = 0; j < numDOF; j += 10) {
                    if (j > 0 && j % (nCell * 10) == 0) {
                        j += nCell * 10;
                        if (j >= numDOF) {
                            break;
                        }
                    }
                    for (int k = 0; k < 10; k++) {
                        newU[j + k] = u[k];
                    }
                }

                for (int j = nCell * 10; j < numDOF; j += 10) {
                    if (j > nCell * 10 && j % (nCell * 10) == 0) {
                        j += nCell * 10;
                        if (j >= numDOF) {
                            break;
                        }
                    }
                    for (int k = 0; k < 10; k++) {
                        newU[j + k] = u[k + 10];
                    }
                }
            }

            coordinateDef.setToU(box.getMoleculeList(), newU);
            coordinateDef.initNominalU(box.getMoleculeList());

        }

        box.setBoundary(boundary);
        double rCScale = 0.475;
        double rC = aDim * nCell * rCScale;
        //System.out.println("Truncation Radius (" + rCScale +" Box Length): " + rC);

        potential = new P2Nitrogen(space, rC);
        potential.setBox(box);

        potentialMaster.addPotential(potential, new ISpecies[]{species, species});

        int xSites = 2 * nCell + 1;
        int ySites = 4 * nCell + 1;
        int zSites = 2 * nCell + 1;
        pairMatrix = new double[xSites][ySites][zSites][4][4][3][3];

        findPair = new FindPairMoleculeIndexBetaN2(space, coordinateDef);
    }
	
	public double[][] get2ndDerivative(int molec0){
	
		int numMolecule = box.getMoleculeList().size();
		int dofTrans = 3;
		double[][] array = new double[dofTrans][dofTrans*numMolecule];
	
		
		DataTensor transTensor = new DataTensor(space);
		MoleculePair pair = new MoleculePair();	
		/*
		 *	Constructing the upper diagonal of the matrix
		 *	(Skipping the molec1 == molec2) 
		 */
		IMolecule molecule0 = coordinateDef.getBox().getMoleculeList().get(molec0);
		pair.mol0 = molecule0;
		
		boolean isReverseOrder = false;
		for(int molec1=0; molec1<numMolecule; molec1++){
			if(molec0 == molec1) continue;
			
			if(molec0 > molec1){
				isReverseOrder = true;
			} else {
				isReverseOrder = false;
			}
			
			
			// Analytical calculation for 3x3 Translational second Derivative
			pair.mol1 = coordinateDef.getBox().getMoleculeList().get(molec1);
		
			int[] index = findPair.getPairMoleculesIndex(pair.mol0, pair.mol1, isReverseOrder);
			boolean isNewPair = findPair.getIsNewPair(index);
			
			if(isNewPair){
				transTensor.E(potential.secondDerivative(pair));
				for(int i=0; i<3; i++){
					for(int j=0; j<3; j++){
						array[i][molec1*dofTrans + j] = transTensor.x.component(i, j);
						pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j] = array[i][molec1*dofTrans + j];
					}
				}
				
				findPair.updateNewMoleculePair(index);
			} else {
				if(isReverseOrder){
					for(int i=0; i<3; i++){
						for(int j=0; j<3; j++){
							array[i][molec1*dofTrans + j] 
							                            = pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][j][i];
						}
					}
				} else {
					for(int i=0; i<3; i++){
						for(int j=0; j<3; j++){
							array[i][molec1*dofTrans + j] 
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
    					array[i][molec0*dofTrans + j] -= array[i][molec1*dofTrans + j] ;         				
	    				
    				}
    			}
    		}
			
			for(int i=0; i<3; i++){
   				for (int j=0; j<3; j++){
   					pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j] = array[i][molec0*dofTrans + j];
   				}
			}
			
    	
		} else {
				
			for(int i=0; i<3; i++){
   				for (int j=0; j<3; j++){
   					array[i][molec0*dofTrans + j]= pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j];
   				}
			}
			
		}
		return array;
		
	}

	public void constructHessianMatrix(int nCell){
		
		int numMolecules = nCell*nCell*nCell*2;
		int interval = nCell*2;
		int dof = 3;
		double[][][] array = new double[interval][dof][numMolecules*dof];
	
			
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
		
		int nCell =2;
		String densityStr = "0.02330";
		
		if(args.length > 0){
			nCell = Integer.parseInt(args[0]);
		}
		
		if(args.length > 1){
			densityStr = args[1];
		}
		int numMolecule = nCell*nCell*nCell*2;
//		System.out.println("Running simulation to construct Hessian Matrix for beta-phase nitrogen");
//		System.out.println("with numMolecule of "+numMolecule + " at density of " + density);
		
		HarmonicBetaNitrogenModelPairMoleculeSequentialHalf2 test = new HarmonicBetaNitrogenModelPairMoleculeSequentialHalf2(Space3D.getInstance(3), numMolecule, densityStr);
		
		long startTime = System.currentTimeMillis();
		String filename = new String ("beta"+numMolecule+"_2ndDer_d"+densityStr+"_new");

		//test.constructHessianMatrix(filename, nCell);
		test.constructHessianMatrix(nCell);
		long endTime = System.currentTimeMillis();
		//System.out.println("Time taken (s): " + (endTime-startTime)/1000);
	}
	
	
	protected Box box;
	protected Space space;
	protected P2Nitrogen potential;
	protected CoordinateDefinitionNitrogen coordinateDef;
	protected PotentialMaster potentialMaster;
	protected double[][][][][][][] pairMatrix;
	protected FindPairMoleculeIndexBetaN2 findPair;
	private static final long serialVersionUID = 1L;
}
