/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.molecule.IMoleculeList;
import etomica.normalmode.ArrayReader1D;
import etomica.normalmode.BasisBigCell;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.Degree;

import java.io.FileWriter;
import java.io.IOException;

/**
 * 
 * This class is written to minimize the lattice energy of the beta-phase
 *  Nitrogen, which has HCP packing with 2 basis atoms using line minimization
 *  
 *  The varying parameters are all the translational d.o.f of the molecules
 *    NO rotational d.o.f.   
 * 
 *
 * @author Tai Boon Tan
 *
 */
public class MinimizeBetaNitrogenTranslationDOF extends Simulation {
	

	public MinimizeBetaNitrogenTranslationDOF(Space space, int[] nC, double density, String fname, double tol){
		super(space);
		this.space = space;
		this.density = density;
		this.nC = nC;
		this.fname = fname;
		this.tolerance = tol;
		
    	energy = new double[3];
    	allValue = new double[3];
		
		double ratio = 1.631;//Math.sqrt(8.0/3.0);
		aDim = Math.pow(4.0/(Math.sqrt(3.0)*ratio*density), 1.0/3.0);
		cDim = aDim*ratio;
		numMolecule = nC[0]*nC[1]*nC[2]*2;
		
		potentialMaster = new PotentialMaster();
		Basis basisHCP = new BasisHcp();
		BasisBigCell basis = new BasisBigCell(space, basisHCP, nC);
		
		species = new SpeciesN2(space);
		addSpecies(species);
		
		box = new Box(space);
		addBox(box);
		box.setNMolecules(species, numMolecule);
		
		boxDim = new Vector[3];
		boxDim[0] = space.makeVector(new double[]{nC[0]*aDim, 0, 0});
		boxDim[1] = space.makeVector(new double[]{-nC[1]*aDim*Math.cos(Degree.UNIT.toSim(60)), nC[1]*aDim*Math.sin(Degree.UNIT.toSim(60)), 0});
		boxDim[2] = space.makeVector(new double[]{0, 0, nC[2]*cDim});
		
		boundary = new BoundaryDeformablePeriodic(space, boxDim);
		primitive = new PrimitiveHexagonal(space, nC[0]*aDim, nC[2]*cDim);
		
		coordinateDefinition = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordinateDefinition.setIsGamma();
		coordinateDefinition.setOrientationVectorGamma(space);
		coordinateDefinition.initializeCoordinates(new int[]{1,1,1});
		
		box.setBoundary(boundary);
		double rC = box.getBoundary().getBoxSize().getX(0)*0.475;
		//System.out.println("rC: " + rC);
		
		P2Nitrogen potential = new P2Nitrogen(space, rC);
		potential.setBox(box);
		
		potentialMaster.addPotential(potential, new ISpecies[]{species, species});
		
		meterPotential = new MeterPotentialEnergy(potentialMaster, box);
	}
	

	public double getEnergy (double[] u){
		
		int numCells =  coordinateDefinition.getBasisCells().length;
		
		for (int cell=0; cell<numCells; cell++){
			IMoleculeList molecules = coordinateDefinition.getBasisCells()[cell].molecules;
			coordinateDefinition.setToU(molecules, u);
		}
    
		return meterPotential.getDataAsScalar();
	}
	

	public void doFindMinimum(double[] parameter){
		this.parameters = parameter;
		int numIter = 1;
		double initEnergy = 0.0;
		double afterEnergy = 0.0;
		double initParam = 0.0;
		double[] minVal = new double[parameter.length];
		double[] maxVal = new double[parameter.length];
	
		for (int i=0; i<parameter.length; i++){
			minVal[i] = parameter[i] *(1 - 0.0011);
			maxVal[i] = parameter[i] *(1 + 0.0012);
		}
		
		while(numIter < 100){
		
			for (int iVar=0; iVar<parameter.length; iVar++){
		//	for(int iVar=parameter.length-1; iVar>0; iVar--){
				if(iVar>0 && (iVar%5==3 || iVar%5==4)){
					continue;
				}
				initEnergy = getEnergy(parameters);
				initParam = parameters[iVar];
				parameters[iVar] = findOptParameter(minVal[iVar], maxVal[iVar], parameters, iVar);
				
				afterEnergy = getEnergy(parameters);
				//System.out.println("diff: "+ (afterEnergy-initEnergy));
				
				if(afterEnergy < initEnergy){
				
		            if(Math.abs(parameters[iVar]) < 1e-8){
		            	minVal[iVar] = parameter[iVar] - 0.01;
						maxVal[iVar] = parameter[iVar] + 0.01;
			        	
		            } else{
						minVal[iVar] = parameter[iVar] *(1-0.003/Math.sqrt(numIter));
						maxVal[iVar] = parameter[iVar] *(1+0.003/Math.sqrt(numIter));
			        
		            }
		 
					//System.out.println("Write FILE DONE");
		            
				} else {
					parameters[iVar] = initParam;
					
					if(Math.abs(afterEnergy - initEnergy)< tolerance){
						System.out.println("Minimum found! with tolerance of " + tolerance);
						System.out.println("after: "+ afterEnergy/numMolecule + " before:" + initEnergy/numMolecule);
						return;
					}
				}
			}
			System.out.println("after " +numIter + " loop; energy: "+ afterEnergy/numMolecule);
			
			++numIter;
			
			try {
				FileWriter fileWriter = new FileWriter(fname+".out", false);
				for (int i=0; i<parameters.length; i++){
					fileWriter.write(parameters[i] + " ");
					
					if(i>1 && i%10==9){
						fileWriter.write("\n");
					}
				}
				
				fileWriter.close();
				
			} catch(IOException e){
				throw new RuntimeException("Failed to write coord data normalize coord U" + e);
			
			}
		}
	
	}
	
	public double findOptParameter(double min, double max, double[] param, int iVar){
    	
		/*
		 * Every time you comes in, you only vary a parameter, called "value"
		 * 
		 * Need to vary min and max
		 */
		
    	int bootstrap = 0;

    	double value = min;
    	param[iVar] = value;
    	
        while (true) {
        	
            latticeEnergy = getEnergy(param);
//          System.out.println("<findOpt> ivar: "+iVar +" ; " +value+" ;lattice energy: " + latticeEnergy/numMolecule);
//			System.out.println(allValue[0]+" "+allValue[1]+" "+allValue[2]);
//			System.out.println(energy[0]+" "+energy[1]+" "+energy[2]);
           
			if (bootstrap < 3) {
                allValue[bootstrap] = value;
                energy[bootstrap] = latticeEnergy;
                bootstrap++;
                
                value += 0.5*(max-min);
                param[iVar] = value;
               
            }
            else {
                if (value > allValue[2]) {
                    allValue[0] = allValue[1];
                    allValue[1] = allValue[2];
                    allValue[2] = value;
                    energy[0] = energy[1];
                    energy[1] = energy[2];
                    energy[2] = latticeEnergy;
                }
                else if (value < allValue[0]) {
                    allValue[2] = allValue[1];
                    allValue[1] = allValue[0];
                    allValue[0] = value;
                    energy[2] = energy[1];
                    energy[1] = energy[0];
                    energy[0] = latticeEnergy;
                }
                else if (energy[2] > energy[0]) {
                    max = allValue[2];
                    if (value > allValue[1]) {
                        energy[2] = latticeEnergy;
                        allValue[2] = value;
                    }
                    else {
                        energy[2] = energy[1];
                        allValue[2] = allValue[1];
                        energy[1] = latticeEnergy;
                        allValue[1] = value;
                    }
                }
                else {
                    min = allValue[0];
                    if (value < allValue[1]) {
                        energy[0] = latticeEnergy;
                        allValue[0] = value;
                    }
                    else {
                        energy[0] = energy[1];
                        allValue[0] = allValue[1];
                        energy[1] = latticeEnergy;
                        allValue[1] = value;
                    }
                }

                if ((energy[1] > energy[0] && energy[1] > energy[2]) || (energy[1] < energy[0] && energy[1] < energy[2])){
                    // we found a maximum, due to numerical precision failure
                    // just bail and pretend that the middle point is the global minimum
                    return allValue[1];
                }
            }

            if (bootstrap == 3) {
                // now estimate minimum in U from the three points.
            	
                double dc01 = allValue[1]-allValue[0];
                double dc12 = allValue[2]-allValue[1];
                double du01 = energy[1]-energy[0];
                double du12 = energy[2]-energy[1];
                double dudc01 = du01/dc01;
                double dudc12 = du12/dc12;
                double m = (dudc12-dudc01)/(0.5*(dc01+dc12));
                value = 0.9*(0.5*(allValue[1]+allValue[2]) - dudc12/m) + 0.1*(0.5*(allValue[0]+allValue[2]));
                param[iVar] = value;
                
                if (value == allValue[1] || value == allValue[2]) {
                    value = 0.5*(allValue[1] + allValue[2]);
                    param[iVar] = value;
                }
                if (value == allValue[0] || value == allValue[1]) {
                    value = 0.5*(allValue[1] + allValue[0]);
                    param[iVar] = value;
                }
                if (value < min) {
                    value = 0.5*(min + allValue[0]);
                    param[iVar] = value;
                }
                if (value > max) {
                    value = 0.5*(max + allValue[2]);
                    param[iVar] = value;
                }
                        
                if (value == allValue[0] || value == allValue[1] || value == allValue[2] ) {
                    // we converged value to numerical precision.
                    //System.out.println("value "+ value);
                	
                	//System.out.println("***"+value + " ;a: " + a+ " ;c: " + c);
                    return value;
                }
                if (Math.abs(energy[0]-energy[1])<1e-10 || Math.abs(energy[1]-energy[2])<1e-10 ||Math.abs(energy[0]-energy[2])<1e-10) {
                	
                	if(energy[0]< energy[1]){
                		value = allValue[0];
                	} else {
                		value = allValue[1];
                	}
                	
                	if(energy[1]< energy[2]){
                		value = allValue[1];
                	} else {
                		value = allValue[2];
                	}
                	
                	if(energy[0]< energy[2]){
                		value = allValue[0];
                	} else {
                		value = allValue[2];
                	}
                	return value;
                }
            }
        }
    }
	
	
	public double getLatticeEnergy(){
		return latticeEnergy;
	}
	
	public static void main(String[] args){
		double density = 0.025;
		int nCell = 6;
		double tol = 1e-10;
		
		if(args.length > 0){
			nCell = Integer.parseInt(args[0]);
		}
		if(args.length > 1){
			tol = Double.parseDouble(args[1]);
		}
		
		
		int[] nC = new int[]{nCell, nCell, nCell};
		String fname = "parameterN2beta_nCell"+nC[0];
		
		double[][] vals = ArrayReader1D.getFromFile(fname+".in");
		double[] parameters = new double[vals.length*vals[0].length];
		
		for(int i=0; i<vals.length; i++){
			for (int j=0; j<vals[0].length; j++){
				parameters[i*vals[0].length + j] = vals[i][j];
			}
		}
		
		System.out.println("Running lattice energy minimization algorithm for beta-N2");
		System.out.println("with nCell: " + nCell + " in each dimension at density of " + density);
		System.out.println("with tolerance of " + tol);
		System.out.println("output file to: " + fname + ".out\n");
		
		MinimizeBetaNitrogenTranslationDOF func = new MinimizeBetaNitrogenTranslationDOF(Space.getInstance(3), nC, density, fname, tol);
		System.out.println("Initial Minimum lattice energy (per molecule): " + func.getEnergy(parameters)/func.numMolecule);
		
		func.doFindMinimum(parameters);
		System.out.println("Minimum lattice energy (per molecule): " + func.getEnergy(func.parameters)/func.numMolecule);
		
		try {
			FileWriter fileWriter = new FileWriter(fname+".out", false);
			for (int i=0; i<parameters.length; i++){
				fileWriter.write(parameters[i] + " ");
				
				if(i>1 && i%10==9){
					fileWriter.write("\n");
				}
			}
			
			fileWriter.close();
			
		} catch(IOException e){
			throw new RuntimeException("Failed to write coord data normalize coord U" + e);
		
		}
	}
	
	
	protected CoordinateDefinitionNitrogen coordinateDefinition;
	protected MeterPotentialEnergy meterPotential;
	protected PotentialMaster potentialMaster;
	protected Box box;
	protected SpeciesN2 species;
	protected double density;
	protected Space space;
	protected Primitive primitive;
	protected Boundary boundary;
	protected Vector[] boxDim;
	protected double aDim, cDim;
	protected double [] parameters;  
	protected int[] nC;
	protected String fname;
 	protected double[] energy;
	protected double[] allValue;
	protected double tolerance;
	
	protected double latticeEnergy;
	protected int numMolecule;
	private static final long serialVersionUID = 1L;
}
