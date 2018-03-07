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
import etomica.normalmode.BasisBigCell;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.Degree;

/**
 * 
 *  This class is written to minimize the lattice energy of the beta-phase
 *  Nitrogen, which has HCP packing with 2 basis atoms using line minimization
 *  
 *  The varying parameters are: 
 *  unit cell A: molecule 1 (u[0] through u[4]) and molecule 2 (u[5] through u[9])
 *  unit cell B: molecule 1 (u[0] through u[4]) and molecule 2 (u[5] through u[9])
 *   - 3 translational + 2 rotational
 *  There are total of 20 degrees of freedom.   
 * 
 *
 * @author Tai Boon Tan
 *
 */
public class MinimizeBetaNitrogenLatticeParameter extends Simulation {
	

	public MinimizeBetaNitrogenLatticeParameter(Space space, int[] nC, double density) {
		super(space);
		this.space = space;
		this.density = density;
		this.nC = nC;

		double ratio = 1.631;
		aDim = Math.pow(4.0 / (Math.sqrt(3.0) * ratio * density), 1.0 / 3.0);
		cDim = aDim * ratio;
		numMolecule = nC[0] * nC[1] * nC[2] * 2;

        species = new SpeciesN2(space);
        addSpecies(species);

		potentialMaster = new PotentialMaster();
		Basis basisHCP = new BasisHcp();
		BasisBigCell basis = new BasisBigCell(space, basisHCP, nC);

		boxDim = new Vector[3];
		boxDim[0] = space.makeVector(new double[]{nC[0] * aDim, 0, 0});
		boxDim[1] = space.makeVector(new double[]{-nC[1] * aDim * Math.cos(Degree.UNIT.toSim(60)), nC[1] * aDim * Math.sin(Degree.UNIT.toSim(60)), 0});
		boxDim[2] = space.makeVector(new double[]{0, 0, nC[2] * cDim});
		boundary = new BoundaryDeformablePeriodic(space, boxDim);
		box = this.makeBox(boundary);
		box.setNMolecules(species, numMolecule);

		primitive = new PrimitiveHexagonal(space, nC[0] * aDim, nC[2] * cDim);

		coordinateDefinition = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordinateDefinition.setIsGamma();
		coordinateDefinition.setOrientationVectorGamma(space);
		coordinateDefinition.initializeCoordinates(new int[]{1, 1, 1});
		double rC = box.getBoundary().getBoxSize().getX(0) * 0.475;
		//System.out.println("rC: " + rC);

		P2Nitrogen potential = new P2Nitrogen(space, rC);
		potential.setBox(box);

		potentialMaster.addPotential(potential, new ISpecies[]{species, species});

		meterPotential = new MeterPotentialEnergy(potentialMaster, box);
	}
	

	public double getEnergy (double[] u){
		
		int numCells =  coordinateDefinition.getBasisCells().length;
		int numDOF = coordinateDefinition.getCoordinateDim();
		int dofPerMol = numDOF/numMolecule; 
		int nCelldofinZ = nC[2]*2*dofPerMol;
		double[] newU = new double[numDOF];

		boolean isUnitCellA = true;
		
		if(true){
			int counter = 1;
			for(int i=0; i<newU.length; i+=10){
				
				if(i>0 && i%nCelldofinZ == 0){
					isUnitCellA = !isUnitCellA;
					
					if((nC[1]%2 == 1) && (i==nCelldofinZ*nC[1]*counter)){
						isUnitCellA = !isUnitCellA;
						++ counter;
					}
				}
				
				if(isUnitCellA){
					newU[i] = u[1];
					newU[i+1] = u[2];
					newU[i+2] = u[3];
					newU[i+3] = u[4];
					newU[i+4] = u[5];
				
					newU[i+5] = u[6];
					newU[i+6] = u[7];
					newU[i+7] = u[8];
					newU[i+8] = u[9];
					newU[i+9] = u[10];
				
				} else {
					newU[i] = u[11];
					newU[i+1] = u[12];
					newU[i+2] = u[13];
					newU[i+3] = u[14];
					newU[i+4] = u[15];
				
					newU[i+5] = u[16];
					newU[i+6] = u[17];
					newU[i+7] = u[18];
					newU[i+8] = u[19];
					newU[i+9] = u[20];
				}
			}
		}


		for (int cell=0; cell<numCells; cell++){
			IMoleculeList molecules = coordinateDefinition.getBasisCells()[cell].molecules;
			coordinateDefinition.setToU(molecules, newU);
		}
    
		
		return meterPotential.getDataAsScalar();
	}
	

	public void doFindMinimum(double[] minVal, double[] maxVal, double[] parameter){
		this.parameters = parameter;
		int numIter = 1;
		double initEnergy = 0.0;
		double afterEnergy = 0.0;
		double initParam = 0.0;
		
		while(numIter < 100){
		
			for (int iVar=1; iVar<parameter.length; iVar++){
		//	for(int iVar=parameter.length-1; iVar>0; iVar--){
				initEnergy = getEnergy(parameters);
				initParam = parameters[iVar];
				
//				System.out.println(parameters[0]+", "+parameters[1]+", "+parameters[2]+", "+parameters[3]+", "+parameters[4]+", "
//				              +"\n "+parameters[5]+", "+parameters[6]+", "+parameters[7]+", "+parameters[8]);
				
				System.out.println(parameters[0]+", "+parameters[1]+", "+parameters[2]+", "+parameters[3]+", "+parameters[4]+", "+parameters[5]+", "
	              +"\n "+parameters[6]+", "+parameters[7]+", "+parameters[8]+", "+parameters[9]+", "+parameters[10]+", "
	              +"\n "+parameters[11]+", "+parameters[12]+", "+parameters[13]+", "+parameters[14]+", "+parameters[15]+", "
	              +"\n "+parameters[16]+", "+parameters[17]+", "+parameters[18]+", "+parameters[19]+", "+parameters[20]+", ");
				System.out.println("initEnergy: " + initEnergy/numMolecule);
				parameters[iVar] = findOptParameter(minVal[iVar], maxVal[iVar], parameters, iVar);
				
				afterEnergy = getEnergy(parameters);
				System.out.println("afterEnergy: " + afterEnergy/numMolecule);
					
				
				
				if(afterEnergy < initEnergy){
				
		            if(Math.abs(parameters[iVar]) < 1e-8){
		            	minVal[iVar] = parameter[iVar] - 0.1;
						maxVal[iVar] = parameter[iVar] + 0.1;
			        	
		            } else{
						minVal[iVar] = parameter[iVar] *(1-0.5/Math.sqrt(numIter));
						maxVal[iVar] = parameter[iVar] *(1+0.5/Math.sqrt(numIter));
			        
		            }
				} else {
					parameters[iVar] = initParam;
					
				}
			}
			
			++numIter;
		}
		System.out.println(parameters[0]+" "+parameters[1]+" "+parameters[2]+" "+parameters[3]+" "+parameters[4]);
			
		//System.out.println("value: " + value + " ; latticeEnergy: " + latticeEnergy/numMolecule);
		
	}
	
	public double findOptParameter(double min, double max, double[] param, int iVar){
    	
		/*
		 * Every time you comes in, you only vary a parameter, called "value"
		 * 
		 * Need to vary min and max
		 */
		
    	int bootstrap = 0;
    	energy = new double[3];
    	allValue = new double[3];
       	
    	double value = min;
    	param[iVar] = value;
    	
        while (true) {
        	
            latticeEnergy = getEnergy(param);
            System.out.println("<findOpt> ivar: "+iVar +" ; " +value+" ;lattice energy: " + latticeEnergy/numMolecule);
			
           // System.exit(1);
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
                        
                if (value == allValue[0] || value == allValue[1] || value == allValue[2]) {
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
		int[] nC = new int[]{6,6,6};
		double density = 0.025;
		/*
		 * with 10 parameters
		 */
//		double[] valMin = new double[]{1.631,  0.525, -0.358,  -0.76447, 1.1032, 
//				                              -0.76461, -1.1610,   0.537,  -0.022};
//		double[] valMax = new double[]{1.631,  0.636, -0.351, -0.76441,  1.251,
//				 						      -0.76441, -1.1602,  0.595,  0.097};
//		double[] parameters = new double[]{1.631, 0.5739673932075715, -0.03548346489501618, -0.7644251420509236, 1.1607890301222756, 
//				 -0.764479724569394, -1.1606699722156573, 0.574052007093417, 0.0357423850524004};
		/*
		 * nA=432
		 * 0.5 * box length
		 */
		double[] valMin = new double[]{1.631, -0.0179, -0.044, -0.0001,  0.565, -0.0515,  
				0.017,  0.026, -0.0001, -0.796, 1.146, 
                0.017, -0.042, -0.0004, -0.796, -1.153,   
               -0.0179,  0.028, -0.0002, 0.565,   0.0495};
		double[] valMax = new double[]{1.631, -0.0173, -0.028, 0.00001, 0.595, -0.0495, 
				0.019,  0.042, 0.0001, -0.741,  1.153,
			    0.019, -0.026, 0.0001, -0.741, -1.146,  
			   -0.0173,  0.044,  0.0001,  0.595,  0.0515};
//		double[] parameters = new double[]{1.631, -0.01743687452578401, -0.03773685126445043, 8.504149130176561E-6, 0.576433699983271, -0.05101361014034025, 
//				 0.01886770973548305, 0.037963868640833384, 2.7547951565875963E-5, -0.7619445259390548, 1.149734539419814, 
//				 0.018503328245454753, -0.037739156902926424, -1.43806022202562E-5, -0.7620319487586951, -1.1496948988796698, 
//				 -0.01778226339270643, 0.03796606350686372, 4.786431025978644E-6, 0.5764589036086002, 0.05087509389025041};

		/*
		 * nA=432
		 * 0.475 * box length
		 */
//		double[] valMin = new double[]{1.631, -0.0379, -0.074, -0.0001,  0.365, -0.1515,  
//				0.007, -0.026, -0.0001, -0.996, 0.846, 
//                0.007, -0.042, -0.0004, -0.996, -1.353,   
//               -0.0379, -0.028, -0.0002, 0.365,   0.1495};
//		double[] valMax = new double[]{1.631, -0.0073, -0.008, 0.00001, 0.795, -0.1495, 
//				0.039,  0.142, 0.0001, -0.341,  1.353,
//			    0.039, -0.026, 0.0001, -0.341, -0.846,  
//			   -0.0073,  0.144,  0.0001,  0.795,  0.1515};
		double[] parameters = new double[]{1.631, -0.017515146245074477, -0.03773538976378145, 6.287613884003924E-6, 0.5764933481704314, -0.05097996592367698, 
				 0.018790453472979852, 0.03796913106646624, 2.5357492230437167E-5, -0.7620276221224523, 1.1497081557927014, 
				 0.018425506729840176, -0.037737747877783084, -1.705003526285108E-5, -0.7621162163806268, -1.1496677505985178, 
				 -0.01785975369940704, 0.037971441969044534, 2.0188668491418435E-6, 0.5765191980050342, 0.05084130993590332};

		
		MinimizeBetaNitrogenLatticeParameter func = new MinimizeBetaNitrogenLatticeParameter(Space.getInstance(3), nC, density);
		
		func.doFindMinimum(valMin, valMax, parameters);
		
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
	protected double[] energy;
	protected double[] allValue;
	
	protected double latticeEnergy;
	protected int numMolecule;
	private static final long serialVersionUID = 1L;
}
