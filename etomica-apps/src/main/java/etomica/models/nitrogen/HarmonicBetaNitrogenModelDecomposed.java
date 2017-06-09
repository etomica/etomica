/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import etomica.box.Box;
import etomica.data.types.DataTensor;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.molecule.MoleculePair;
import etomica.normalmode.ArrayReader1D;
import etomica.normalmode.BasisBigCell;
import etomica.normalmode.CoordinateDefinition;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.units.Degree;

import java.io.FileWriter;
import java.io.IOException;



/**
 * 
 * Class that constructs the Hessian Matrix
 * 
 * Special features about this class is that: 
 * 	1. each row of blocks can be calculated in separate simulations 
 * 	   (to account for the expensive nature of constructing the Matrix in a single run)
 *  2. you can combine all the row blocks into the (N*dof) x (N*dof) matrix once all the simulations are done
 *     N is the number of molecules and dof is the degrees of freedom for each molecule
 * 
 * @author Tai Boon Tan
 *
 */
public class HarmonicBetaNitrogenModelDecomposed extends Simulation{

	
	public HarmonicBetaNitrogenModelDecomposed(Space space, int numMolecule, double density) {
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
		
		coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordinateDef.setIsBeta();
		coordinateDef.setOrientationVectorBeta(space);
		coordinateDef.initializeCoordinates(nCells);
		
		box.setBoundary(boundary);
		double rCScale = 0.475;
		double rC =aDim*nC*rCScale;
		System.out.println("Truncation Radius (" + rCScale +" Box Length): " + rC);
		
		potential = new P2Nitrogen(space, rC);
		potential.setBox(box);

		potentialMaster.addPotential(potential, new ISpecies[]{species, species});

	}
	
	
	
	public double[][] get2ndDerivativeMoleculei(CoordinateDefinition coordinateDef, int iMolecule){
		
		int numMolecule = box.getMoleculeList().getMoleculeCount();
		int dofPerMol = coordinateDef.getCoordinateDim()/numMolecule;
		
		double[][] array = new double[dofPerMol][coordinateDef.getCoordinateDim() - ((iMolecule+1)*dofPerMol)];
		double[] newU = new double[coordinateDef.getCoordinateDim()];
		DataTensor transTensor = new DataTensor(space);
		MoleculePair pair = new MoleculePair();
		
		CalcNumerical2ndDerivative cm2ndD = new CalcNumerical2ndDerivative(box, potentialMaster, coordinateDef);
		
		/*
		 *	Constructing the upper diagonal of the matrix
		 *	(Skipping the molec1 == molec2) 
		 */
	
		pair.atom0 = box.getMoleculeList().getMolecule(iMolecule);
			
		for(int molec1=(iMolecule+1); molec1<numMolecule; molec1++){
			/*
			 * working within the 5x5 Matrix
			 */
			// Analytical calculation for 3x3 Translational second Derivative
			pair.atom1 = box.getMoleculeList().getMolecule(molec1);
		
			transTensor.E(potential.secondDerivative(pair));
			for(int i=0; i<3; i++){
				for(int j=0; j<3; j++){
					array[i][(molec1-(iMolecule+1))*dofPerMol + j] = transTensor.x.component(i, j);
				}
			}
			// Numerical calculation for the Cross (trans and rotation) and rotation second Derivative
			for(int i=0; i<dofPerMol; i++){
				for(int j=0; j<dofPerMol; j++){
					if(i<3 && j<3) continue;
					array[i][(molec1-(iMolecule+1))*dofPerMol + j] = cm2ndD.d2phi_du2(new int[]{(iMolecule*dofPerMol + i),(molec1*dofPerMol + j)}, newU);
			
				}
			}
		}
		return array;
		
	}
	
	public double[][] contructFullMatrix(CoordinateDefinition coordinateDef, String fname){
		
		double[][] array = new double[coordinateDef.getCoordinateDim()][coordinateDef.getCoordinateDim()];
		double[] newU = new double[coordinateDef.getCoordinateDim()];
		
		int numMolecule = box.getMoleculeList().getMoleculeCount();
		int dofPerMol = coordinateDef.getCoordinateDim()/numMolecule;
		
		CalcNumerical2ndDerivative cm2ndD = new CalcNumerical2ndDerivative(box, potentialMaster, coordinateDef);
		
		/*
		 *	Constructing the upper diagonal of the matrix
		 *	(Skipping the molec1 == molec2) 
		 */
		for(int molec0=0; molec0<(numMolecule-1); molec0++){
			double[][] tempArray = getMatrixFromFile(fname, molec0);
			
			for(int molec1=(molec0+1); molec1<numMolecule; molec1++){
				for(int i=0; i<dofPerMol; i++){
					for (int j=0; j<dofPerMol; j++){
						array[molec0*dofPerMol + i][molec1*dofPerMol + j] = tempArray[i][(molec1-(molec0+1))*dofPerMol + j];
						
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
			for(int molec1=0; molec1<numMolecule; molec1++){
    			if(molec0==molec1) continue; // we might double sum the elements in array[a][a] if we don't skip the pair
    			for(int i=0; i<3; i++){
    				for (int j=0; j<3; j++){
    					array[molec0*dofPerMol + i][molec0*dofPerMol + j] -= array[molec0*dofPerMol + i][molec1*dofPerMol + j] ;         				
    				}
    			}
    		}
    	}
		
		for(int molec0=0; molec0<numMolecule; molec0++){
			// Numerical calculation for the Cross (trans and rotation) and rotation second Derivative
			for(int i=0; i<dofPerMol; i++){
				for(int j=0; j<dofPerMol; j++){
					if(i<3 && j<3) continue;
					array[molec0*dofPerMol + i][molec0*dofPerMol + j] = cm2ndD.d2phi_du2(new int[]{(molec0*dofPerMol + i),(molec0*dofPerMol + j)}, newU);
				}    		
    		}
    	}
		
		return array;
		
	}
	
	public double[][] getMatrixFromFile(String fname, int iMolecule){
		double[][] array = ArrayReader1D.getFromFile(fname+"_mol"+iMolecule);
		return array;
	}
	
	public void doEigenDecompose(double[][] array, String filename){
		
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
	
	public static void main (String[] args){
		
		int numMolecule =432;
		int iMolecule = 430;
		double density = 0.025;
		boolean isCombineFile = false;
		
		if(args.length > 0){
			numMolecule = Integer.parseInt(args[0]);
		}
		if(args.length > 1){
			iMolecule = Integer.parseInt(args[1]);
		}
		if(args.length > 2){
			isCombineFile = Boolean.parseBoolean(args[2]);
		}
		String filename = "beta"+numMolecule+"_2ndDer_d"+density;
		
		System.out.println("Running Hessian Matrix Construction Program for Beta-phase Nitrogen Model");
		System.out.println("with density of " + density);
		System.out.println("isCombineFile: " +  isCombineFile);
		System.out.println("output file to: " + filename);
		HarmonicBetaNitrogenModelDecomposed test = new HarmonicBetaNitrogenModelDecomposed(Space3D.getInstance(3), numMolecule, density);
	
		long startTime = System.currentTimeMillis();
		
		if(isCombineFile){
			double[][] array = test.contructFullMatrix(test.coordinateDef, filename);
			try {
				FileWriter fileWriter = new FileWriter(filename+"_all");
				
				for (int i=0; i<array.length; i++){
					for (int j=0; j<array[0].length; j++){
						double value = array[i][j];
						fileWriter.write(value+ " ");
					}
					fileWriter.write("\n");
				}
				fileWriter.close();
				
			} catch (IOException e) {
				
			}
			
			test.doEigenDecompose(array, filename);
			
		} else {
			System.out.println("Contructing matrix for molecule " + iMolecule);
			double[][]array = test.get2ndDerivativeMoleculei(test.coordinateDef, iMolecule);
			
			String fname = new String (filename+"_mol"+iMolecule);
			try {
				FileWriter fileWriter = new FileWriter(fname);
				
				for (int i=0; i<array.length; i++){
					for (int j=0; j<array[0].length; j++){
						double value = array[i][j];
						fileWriter.write(value+ " ");
					}
					fileWriter.write("\n");
				}
				fileWriter.close();
				
			} catch (IOException e) {
				
			}
		}

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
