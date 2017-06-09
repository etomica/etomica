/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import etomica.api.ISpecies;
import etomica.box.Box;
import etomica.data.types.DataTensor;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.molecule.MoleculePair;
import etomica.normalmode.BasisBigCell;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.units.Degree;

import java.io.FileWriter;
import java.io.IOException;

/**
 * 
 * Beta-phase Nitrogen Harmonic Approximation
 * with only 3 translational degrees of freedom
 * 
 * @author Tai Boon Tan
 *
 */
public class HarmonicBetaNitrogenModel extends Simulation{
	
	public HarmonicBetaNitrogenModel(Space space, int numMolecule, double density) {
		super(space);
		this.space = space;
		
		potentialMaster = new PotentialMaster();
		
	  	double ratio = 1.631;
		double aDim = Math.pow(4.0/(Math.sqrt(3.0)*ratio*density), 1.0/3.0);
		double cDim = aDim*ratio;
		System.out.println("\naDim: " + aDim + " ;cDim: " + cDim);
		int nC = (int)Math.pow(numMolecule/1.999999999, 1.0/3.0);
		
		Basis basisHCP = new BasisHcp();
		BasisBigCell basis = new BasisBigCell(space, basisHCP, new int[]{nC,nC,nC});
        
		Vector[] boxDim = new Vector[3];
		boxDim[0] = space.makeVector(new double[]{nC*aDim, 0, 0});
		boxDim[1] = space.makeVector(new double[]{-nC*aDim*Math.cos(Degree.UNIT.toSim(60)), nC*aDim*Math.sin(Degree.UNIT.toSim(60)), 0});
		boxDim[2] = space.makeVector(new double[]{0, 0, nC*cDim});
		
		int[] nCells = new int[]{1,1,1};
		Boundary boundary = new BoundaryDeformablePeriodic(space, boxDim);
		Primitive primitive = new PrimitiveHexagonal(space, nC*aDim, nC*cDim);
		
		SpeciesN2 species = new SpeciesN2(space);
		addSpecies(species);
		
		box = new Box(space);
		addBox(box);
		box.setNMolecules(species, numMolecule);		
		
		coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space, 0);
		coordinateDef.setIsBeta();
		coordinateDef.setOrientationVectorBeta(space);
		coordinateDef.initializeCoordinates(nCells);
		
		box.setBoundary(boundary);
		double rCScale = 0.475;
		double rC = aDim*nC*rCScale;
		System.out.println("Truncation Radius (" + rCScale +" Box Length): " + rC);
		
		potential = new P2Nitrogen(space, rC);
		potential.setBox(box);

		potentialMaster.addPotential(potential, new ISpecies[]{species, species});

	}
	
	public double[][] get2ndDerivative(){
	
		int numMolecule = box.getMoleculeList().getMoleculeCount();
		int dofTrans = 3;
		double[][] array = new double[dofTrans*numMolecule][dofTrans*numMolecule];
			
		DataTensor transTensor = new DataTensor(space);
		MoleculePair pair = new MoleculePair();	
		/*
		 *	Constructing the upper diagonal of the matrix
		 *	(Skipping the molec1 == molec2) 
		 */
		for(int molec0=0; molec0<numMolecule; molec0++){
			pair.atom0 = box.getMoleculeList().getMolecule(molec0);
			
			for(int molec1=molec0; molec1<numMolecule; molec1++){
				if(molec0 == molec1) continue;
				
				// Analytical calculation for 3x3 Translational second Derivative
				pair.atom1 = box.getMoleculeList().getMolecule(molec1);
		
				transTensor.E(potential.secondDerivative(pair));
				for(int i=0; i<3; i++){
					for(int j=0; j<3; j++){
						array[molec0*dofTrans + i][molec1*dofTrans + j] = transTensor.x.component(i, j);
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
				for(int i=0; i<3; i++){
					for (int j=0; j<3; j++){
						array[molec0*dofTrans + i][molec1*dofTrans + j] = array[molec1*dofTrans + j][molec0*dofTrans + i];
						
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
			for(int molec1=0; molec1<numMolecule; molec1++){
    			if(molec0==molec1) continue; // we might double sum the elements in array[a][a] if we don't skip the pair
    			for(int i=0; i<3; i++){
    				for (int j=0; j<3; j++){
    					array[molec0*dofTrans + i][molec0*dofTrans + j] -= array[molec0*dofTrans + i][molec1*dofTrans + j] ;         				
    				}
    			}
    		}
    	}
		
		return array;
		
	}
	
	public void doEigenDecomposeAndFile(double[][] array, String filename){
		
		try{
			
			FileWriter fileWriterVal = new FileWriter(filename+".val");
			FileWriter fileWriterVec = new FileWriter(filename+".vec");
			
				
			Matrix matrix = new Matrix(array);
			EigenvalueDecomposition ed = matrix.eig();
			double[] eVals = ed.getRealEigenvalues();
			double[][] eVecs = ed.getV().getArray();
				
			// output .val file
			for (int ival=0; ival<eVals.length; ival++){
				if (eVals[ival] < 1E-10){
					fileWriterVal.write("0.0 ");
				} else {
					fileWriterVal.write(1/eVals[ival]+ " ");
                }
			}
			fileWriterVal.write("\n");

			// output .vec file
			for (int ivec=0; ivec<eVecs.length; ivec++ ){
				for(int jvec=0; jvec<eVecs[0].length; jvec++){
					if (Math.abs(eVecs[jvec][ivec])<1e-10){
						fileWriterVec.write("0.0 ");
					} else {
                		fileWriterVec.write(eVecs[jvec][ivec] + " ");
               		}
               	}
               	fileWriterVec.write("\n");
			}
		
			fileWriterVal.close();
            fileWriterVec.close();
		
		} catch (IOException e) {
			
		}
	}
	
	public void get2ndDerivativeFile(double[][] array, String filename){
		
		int dofTrans = array.length;
		try {
			FileWriter fileWriter = new FileWriter(filename);//+".hes");
			
			for (int i=0; i<dofTrans; i++){
				for (int j=0; j<dofTrans; j++){
					double value = array[i][j];
					fileWriter.write(value+ " ");
				}
				fileWriter.write("\n");
			}
			fileWriter.close();
			
		} catch (IOException e) {
			
		}
	}
	
	public static void main (String[] args){
		
		int nCell = 8;
		double density = 0.025;
		
		if(args.length > 0){
			nCell = Integer.parseInt(args[0]);
		}
		int numMolecule = nCell*nCell*nCell*2;
		System.out.println("Running simulation to construct Hessian Matrix for beta-phase nitrogen");
		System.out.println("with numMolecule of "+numMolecule + " at density of " + density);
		
		HarmonicBetaNitrogenModel test = new HarmonicBetaNitrogenModel(Space3D.getInstance(3), numMolecule, density);

		long startTime = System.currentTimeMillis();
		double[][]testArray = test.get2ndDerivative();
		String filename = new String ("beta"+numMolecule+"_2ndDer_d"+density+"_new");

		test.get2ndDerivativeFile(testArray, filename);
		System.out.println("***get second derivative file done!");
//		test.doEigenDecomposeAndFile(testArray, filename);
//		System.out.println("***do decomposition done!");
		
		long endTime = System.currentTimeMillis();
		System.out.println("Time taken (s): " + (endTime-startTime)/1000);
	}
	
	
	protected Box box;
	protected Space space;
	protected P2Nitrogen potential;
	protected CoordinateDefinitionNitrogen coordinateDef;
	protected PotentialMaster potentialMaster;
	private static final long serialVersionUID = 1L;
}
