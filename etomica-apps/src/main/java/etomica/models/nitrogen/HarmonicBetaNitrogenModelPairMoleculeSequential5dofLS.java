/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeChildAtomAction;
import etomica.box.Box;
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
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.units.Degree;

/**
 * 
 * Beta-phase Nitrogen Harmonic Approximation
 * with only 5 degrees of freedom
 * 
 * with lattice sum
 * 
 * @author Tai Boon Tan
 *
 */
public class HarmonicBetaNitrogenModelPairMoleculeSequential5dofLS extends Simulation{
	
	public HarmonicBetaNitrogenModelPairMoleculeSequential5dofLS(Space space, int numMolecule, double density, double rC) {
        super(space);
        this.space = space;

        potentialMaster = new PotentialMaster();

        double ratio = 1.631;
        double aDim = Math.pow(4.0 / (Math.sqrt(3.0) * ratio * density), 1.0 / 3.0);
        double cDim = aDim * ratio;
        //System.out.println("\naDim: " + aDim + " ;cDim: " + cDim);
        int nCell = (int) Math.pow(numMolecule / 1.999999999, 1.0 / 3.0);

        Basis basisHCP = new BasisHcp();
        BasisBigCell basis = new BasisBigCell(space, basisHCP, new int[]{nCell, nCell, nCell});

        Vector[] boxDim = new Vector[3];
        boxDim[0] = space.makeVector(new double[]{nCell * aDim, 0, 0});
        boxDim[1] = space.makeVector(new double[]{-nCell * aDim * Math.cos(Degree.UNIT.toSim(60)), nCell * aDim * Math.sin(Degree.UNIT.toSim(60)), 0});
        boxDim[2] = space.makeVector(new double[]{0, 0, nCell * cDim});

        int[] nCells = new int[]{1, 1, 1};
        Boundary boundary = new BoundaryDeformablePeriodicSwitch(space, boxDim);
        Primitive primitive = new PrimitiveHexagonal(space, nCell * aDim, nCell * cDim);

        SpeciesN2 species = new SpeciesN2(space);
        addSpecies(species);

        box = this.makeBox(boundary);
        box.setNMolecules(species, numMolecule);

        coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
        coordinateDef.setIsBeta();
        coordinateDef.setOrientationVectorBeta(space);
        coordinateDef.initializeCoordinates(nCells);

        double[] u = new double[20];
        if (true) {
            BetaPhaseLatticeParameterLS parameters = new BetaPhaseLatticeParameterLS();
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

        this.rC = rC;
        //System.out.println("Truncation Radius (" + rCScale +" Box Length): " + rC);

        potential = new P2Nitrogen(space, rC);
        potential.setEnablePBC(false);
        potential.setBox(box);

        potentialMaster.addPotential(potential, new ISpecies[]{species, species});

        int xSites = 2 * nCell + 1;
        int ySites = 4 * nCell + 1;
        int zSites = 2 * nCell + 1;
        pairMatrix = new double[xSites][ySites][zSites][4][4][5][5];
        isFoundReverse = new boolean[xSites][ySites][zSites][4][4]; //default to false

        cAN2nD = new CalcHalfAnalyticHalfNumeric2ndDerivativeNitrogen(space, box, potential, coordinateDef, true, rC, false);
        cA2nD = new CalcAnalytical2ndDerivativeNitrogen(space, box, potential, coordinateDef, true, rC);


        findPair = new FindPairMoleculeIndexBetaN2(space, coordinateDef);

        translateBy = new AtomActionTranslateBy(coordinateDef.getPrimitive().getSpace());
        atomGroupActionTranslate = new MoleculeChildAtomAction(translateBy);
        lsPosition = coordinateDef.getPrimitive().getSpace().makeVector();

        xVecBox = Math.sqrt(box.getBoundary().getEdgeVector(0).squared());
        yVecBox = Math.sqrt(box.getBoundary().getEdgeVector(1).squared());
        zVecBox = Math.sqrt(box.getBoundary().getEdgeVector(2).squared());

//		System.out.println(xVecBox+" "+yVecBox+ " " + zVecBox);
        double rX = xVecBox;
        this.nLayer = (int) Math.round(rC / rX + 0.5);

//		System.out.println("rX: " + rX);
//		System.out.println("nLayer: " + nLayer);
//		System.exit(1);
    }
	
	public double[][] get2ndDerivative(int molec0){
	
		int numMolecule = box.getMoleculeList().size();
		int dofPerMol = coordinateDef.getCoordinateDim()/numMolecule;
		double[][] array = new double[dofPerMol][coordinateDef.getCoordinateDim()];
		
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
			
			
			IMolecule molecule1 = coordinateDef.getBox().getMoleculeList().get(molec1);
			pair.mol1 = molecule1;
		
			int[] index = findPair.getPairMoleculesIndex(pair.mol0, pair.mol1, isReverseOrder);
			boolean isNewPair = findPair.getIsNewPair(index);
			
			if(isReverseOrder && isNewPair){
				isFoundReverse[index[0]][index[1]][index[2]][index[3]][index[4]] = true;
			}
			
			if(isNewPair){
				
//				transTensor.E(0.0);
				((BoundaryDeformablePeriodicSwitch)box.getBoundary()).setDoPBC(false);
				
				//do Lattice sum
				double[][] sumA = new double[5][5];
				for(int x=-nLayer; x<=nLayer; x++){
					for(int y=-nLayer; y<=nLayer; y++){
						for(int z=-nLayer; z<=nLayer; z++){
							lsPosition.E(new double[]{x*xVecBox-y*yVecBox*Math.cos(Degree.UNIT.toSim(60)), 
									 y*yVecBox*Math.sin(Degree.UNIT.toSim(60)), 
									 z*zVecBox});
							translateBy.setTranslationVector(lsPosition);
							atomGroupActionTranslate.actionPerformed(molecule1);
		
//							transTensor.PE(potential.secondDerivative(pair));
														
							double[][] a = cA2nD.d2phi_du2(new int[]{molec1,molec0});
							for(int i=0; i<dofPerMol; i++){
								for(int j=0; j<dofPerMol; j++){
									sumA[i][j] += a[i][j];
								}
							}
							
							lsPosition.TE(-1);
							translateBy.setTranslationVector(lsPosition);
							atomGroupActionTranslate.actionPerformed(molecule1);
						}	
					}	
				}
				((BoundaryDeformablePeriodicSwitch)box.getBoundary()).setDoPBC(true);
				
				for(int i=0; i<dofPerMol; i++){
					for(int j=0; j<dofPerMol; j++){
						array[i][molec1*dofPerMol + j] = sumA[i][j];
						pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j] = array[i][molec1*dofPerMol + j];
					}
				}
				
				findPair.updateNewMoleculePair(index);
			} else {
				
				if(isFoundReverse[index[0]][index[1]][index[2]][index[3]][index[4]]==true){
					isReverseOrder = !isReverseOrder; 
				}
				
				if(isReverseOrder){
					for(int i=0; i<5; i++){
						for(int j=0; j<5; j++){
							array[i][molec1*dofPerMol + j] 
							                            = pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][j][i];
						}
					}
				} else {
					for(int i=0; i<5; i++){
						for(int j=0; j<5; j++){
							array[i][molec1*dofPerMol + j] 
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
			((BoundaryDeformablePeriodicSwitch)box.getBoundary()).setDoPBC(false);
			for(int i=0; i<dofPerMol; i++){
				for(int j=0; j<dofPerMol; j++){
					if(i<3 && j<3) continue;
					array[i][molec0*dofPerMol + j] = cAN2nD.d2phi_du2(new int[]{molec0,molec0}, new int[]{i,j});
					pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j] = array[i][molec0*dofPerMol + j];
				}    		
	    	}
						
			((BoundaryDeformablePeriodicSwitch)box.getBoundary()).setDoPBC(true);
			findPair.updateNewMoleculePair(index);
				
		} else {
			for(int i=0; i<dofPerMol; i++){
				for(int j=0; j<dofPerMol; j++){
					if(i<3 && j<3) continue;
					array[i][molec0*dofPerMol + j] =  pairMatrix[index[0]][index[1]][index[2]][index[3]][index[4]][i][j];
				}    		
	   		}
			
		}
		
//		if(molec0==3){
//			for(int i=0; i<dofPerMol; i++){
//				for(int j=0; j<dofPerMol; j++){
//					System.out.print(array[i][molec0*dofPerMol + j] + " ");
//				}    		
//				System.out.println();
//	   		}
//			
//			System.exit(1);
//		}
		
//		for(int molec1=0; molec1<numMolecule; molec1++){
//			if(molec0==molec1) continue; 
//			for(int i=0; i<3; i++){
//				for (int j=0; j<3; j++){
//					array[i][molec0*dofPerMol + j] -= array[i][molec1*dofPerMol + j] ;         				
//    				
//				}
//			}
//		}
//		
//		// Numerical calculation for the Cross (trans and rotation) and rotation second Derivative
//		
//		((BoundaryDeformablePeriodicSwitch)box.getBoundary()).setDoPBC(false);
//		for(int i=0; i<dofPerMol; i++){
//			for(int j=0; j<dofPerMol; j++){
//				if(i<3 && j<3) continue;
//			// j i because it got switched molecule A and molecule B
//			//setting the off-diagonal element to zero
//				array[i][molec0*dofPerMol + j] = cAN2nD.d2phi_du2(new int[]{molec0,molec0}, new int[]{i,j});
//			}    		
//    	}
//					
//		((BoundaryDeformablePeriodicSwitch)box.getBoundary()).setDoPBC(true);
		 
		return array;
		
	}
	
	
	public void constructHessianMatrix(int nCell){
		
		int numMolecules = nCell*nCell*nCell*2;
		int interval = nCell*2;
		int dof = coordinateDef.getCoordinateDim()/numMolecules;
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
		
		int nCell =4;
		double density = 0.0230;
		double rC = 50;
		
		if(args.length > 0){
			nCell = Integer.parseInt(args[0]);
		}
		
		if(args.length > 1){
			density = Double.parseDouble(args[1]);
		}
		if(args.length > 2){
			rC = Double.parseDouble(args[2]);
		}
		
		int numMolecule = nCell*nCell*nCell*2;
//		System.out.println("Running simulation to construct Hessian Matrix for beta-phase nitrogen");
//		System.out.println("with numMolecule of "+numMolecule + " at density of " + density);
		
		HarmonicBetaNitrogenModelPairMoleculeSequential5dofLS test = new HarmonicBetaNitrogenModelPairMoleculeSequential5dofLS(Space3D.getInstance(3), numMolecule, density, rC);
		
		long startTime = System.currentTimeMillis();
		String filename = new String ("beta"+numMolecule+"_2ndDer_d"+density+"_new");

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
	protected CalcHalfAnalyticHalfNumeric2ndDerivativeNitrogen cAN2nD;
	protected CalcAnalytical2ndDerivativeNitrogen cA2nD;
	protected FindPairMoleculeIndexBetaN2 findPair;
	protected AtomActionTranslateBy translateBy;
	protected MoleculeChildAtomAction atomGroupActionTranslate;
	protected Vector lsPosition;
	protected double xVecBox, yVecBox, zVecBox, rC;
	protected int nLayer;
	protected boolean[][][][][] isFoundReverse;
	
	private static final long serialVersionUID = 1L;
}
