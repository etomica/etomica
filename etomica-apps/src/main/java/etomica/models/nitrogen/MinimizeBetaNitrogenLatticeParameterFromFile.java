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
import etomica.math.numerical.ArrayReader1D;
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

import java.io.FileWriter;
import java.io.IOException;

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
public class MinimizeBetaNitrogenLatticeParameterFromFile extends Simulation {
	

	public MinimizeBetaNitrogenLatticeParameterFromFile(Space space, int[] nC, double density, double[] u) {
		super(space);
		this.space = space;
		this.density = density;
		this.nC = nC;

		double ratio = 1.631;
		aDim = Math.pow(4.0 / (Math.sqrt(3.0) * ratio * density), 1.0 / 3.0);
		cDim = aDim * ratio;
		numMolecule = nC[0] * nC[1] * nC[2] * 2;

		potentialMaster = new PotentialMaster();
		Basis basisHCP = new BasisHcp();
		BasisBigCell basis = new BasisBigCell(space, basisHCP, nC);

		species = new SpeciesN2(space);
		addSpecies(species);

		boxDim = new Vector[3];
		boxDim[0] = space.makeVector(new double[]{nC[0] * aDim, 0, 0});
		boxDim[1] = space.makeVector(new double[]{-nC[1] * aDim * Math.cos(Degree.UNIT.toSim(60)), nC[1] * aDim * Math.sin(Degree.UNIT.toSim(60)), 0});
		boxDim[2] = space.makeVector(new double[]{0, 0, nC[2] * cDim});
		boundary = new BoundaryDeformablePeriodic(space, boxDim);
		box = this.makeBox(boundary);
		box.setNMolecules(species, numMolecule);

		primitive = new PrimitiveHexagonal(space, nC[0] * aDim, nC[2] * cDim);

		coordinateDefinition = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordinateDefinition.setIsBeta();
		coordinateDefinition.setOrientationVectorBeta(space);
		coordinateDefinition.initializeCoordinates(new int[]{1, 1, 1});

		int numCells = coordinateDefinition.getBasisCells().length;
		int numDOF = coordinateDefinition.getCoordinateDim();
		int dofPerMol = numDOF / numMolecule;
		int nCelldofinZ = nC[2] * 2 * dofPerMol;
		double[] newU = new double[numDOF];

		boolean isUnitCellA = true;

		if (true) {
			int counter = 1;
			for (int i = 0; i < newU.length; i += 10) {

				if (i > 0 && i % nCelldofinZ == 0) {
					isUnitCellA = !isUnitCellA;

					if ((nC[1] % 2 == 1) && (i == nCelldofinZ * nC[1] * counter)) {
						isUnitCellA = !isUnitCellA;
						++counter;
					}
				}

				if (isUnitCellA) {
					newU[i] = u[0];
					newU[i + 1] = u[1];
					newU[i + 2] = u[2];
					newU[i + 3] = u[3];
					newU[i + 4] = u[4];

					newU[i + 5] = u[5];
					newU[i + 6] = u[6];
					newU[i + 7] = u[7];
					newU[i + 8] = u[8];
					newU[i + 9] = u[9];

				} else {
					newU[i] = u[10];
					newU[i + 1] = u[11];
					newU[i + 2] = u[12];
					newU[i + 3] = u[13];
					newU[i + 4] = u[14];

					newU[i + 5] = u[15];
					newU[i + 6] = u[16];
					newU[i + 7] = u[17];
					newU[i + 8] = u[18];
					newU[i + 9] = u[19];
				}
			}
		}


		for (int cell = 0; cell < numCells; cell++) {
			IMoleculeList molecules = coordinateDefinition.getBasisCells()[cell].molecules;
			coordinateDefinition.setToU(molecules, newU);
		}
		double rC = nC[0] * aDim * 0.475;
		//System.out.println("rC: " + rC);
		P2Nitrogen potential = new P2Nitrogen(space, rC);
		potential.setBox(box);

		potentialMaster.addPotential(potential, new ISpecies[]{species, species});

		meterPotential = new MeterPotentialEnergy(potentialMaster, box);
		initLat = meterPotential.getDataAsScalar() / numMolecule;
		//System.out.println("lattice energy: "+ meterPotential.getDataAsScalar()/numMolecule);
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
					newU[i] = u[0];
					newU[i+1] = u[1];
					newU[i+2] = u[2];
					newU[i+3] = u[3];
					newU[i+4] = u[4];
				
					newU[i+5] = u[5];
					newU[i+6] = u[6];
					newU[i+7] = u[7];
					newU[i+8] = u[8];
					newU[i+9] = u[9];
				
				} else {
					newU[i] = u[10];
					newU[i+1] = u[11];
					newU[i+2] = u[12];
					newU[i+3] = u[13];
					newU[i+4] = u[14];
				
					newU[i+5] = u[15];
					newU[i+6] = u[16];
					newU[i+7] = u[17];
					newU[i+8] = u[18];
					newU[i+9] = u[19];
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
		
		while(numIter < 20){
		
			for (int iVar=0; iVar<parameter.length; iVar++){
		//	for(int iVar=parameter.length-1; iVar>0; iVar--){
				initEnergy = getEnergy(parameters);
				initParam = parameters[iVar];
				

//				System.out.println(parameters[0]+",\n "+parameters[1]+", "+parameters[2]+", "+parameters[3]+", "+parameters[4]+", "+parameters[5]+", "
//	              +"\n "+parameters[6]+", "+parameters[7]+", "+parameters[8]+", "+parameters[9]+", "+parameters[10]+", "
//	              +"\n "+parameters[11]+", "+parameters[12]+", "+parameters[13]+", "+parameters[14]+", "+parameters[15]+", "
//	              +"\n "+parameters[16]+", "+parameters[17]+", "+parameters[18]+", "+parameters[19]+", "+parameters[20]+", ");
//				System.out.println("initEnergy: " + initEnergy/numMolecule);
				parameters[iVar] = findOptParameter(minVal[iVar], maxVal[iVar], parameters, iVar);
				
				afterEnergy = getEnergy(parameters);
				//System.out.println("afterEnergy: " + afterEnergy/numMolecule);
					
				if(afterEnergy < initEnergy && afterEnergy<initLat){
					System.out.println("**************** LOWER ENERGY LATTICE STRUCTURE FOUND! *******************");
					for(int i=0; i<parameter.length; i++){
						System.out.print(parameters[i]+", ");
						if((i+1)%5==0){
							System.out.println("");
						}
					}
					System.out.println(numIter + " "+ iVar+" lower lattice energy (sim unit): " + getEnergy(parameters));
					
					
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
	
//		for (int i=0; i<parameters.length;i++){
//			System.out.println(parameters[i]);
//		}
//		System.out.println(getEnergy(parameters)/numMolecule);
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
            //System.out.println("<findOpt> ivar: "+iVar +" ; " +value+" ;lattice energy: " + latticeEnergy/numMolecule);
			
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
		
		String filename = "/tmp/inputd0.02300";
		double density = 0.02300;
		int nCells = 8;
		double scale = 6e-1;
		
	    if(args.length > 0){
			filename = args[0];
		}
		if(args.length > 1){
			density = Double.parseDouble(args[1]);
		}
		if(args.length > 2){
			nCells = Integer.parseInt(args[2]);
		}
		if(args.length > 3){
			scale = Double.parseDouble(args[3]);
		}
        
		int[] nC = new int[]{nCells,nCells,nCells};
//		double[] parameters = new double[21];
//		double[][] paramFromFile = ArrayReader1D.getFromFile(filename);
//		for (int i=0; i<parameters.length;i++){
//			parameters[i] = paramFromFile[i][0];
//		}
		
		double[] parameters = new double[20];
		double[][] paramFromFile = ArrayReader1D.getFromFile(filename);
		int k=0;
		for (int i=0; i<paramFromFile.length;i++){
			for (int j=0; j<paramFromFile[0].length;j++){
				
				parameters[k]=paramFromFile[i][j];
				k++;
			}	
		}

		//rho = 0.02300
//		double[] parameters = new double[]{
//				0.008643921859739718, 0.0026272419155209253, 0.0019564392850298443, -0.005368679722793812, -0.01963128167700018, 
//				-0.007220434499475046, 7.965716892259778E-4, 0.0019506990169544465, 0.005210028505749975, -0.019695510442517328, 
//				-0.007182981317926607, 0.0026351865272155905, 0.001963353387062189, -0.005023510603676701, -0.019828787213698948, 
//				0.008621853435524681, 7.735415935335443E-4, 0.001965354737210446, 0.005299714904072902, -0.019643975672095665
//		};
		
		//rho = 0.02320
//		double[] parameters = new double[]{
//				0.007559425210441226, 0.002144669727082794, 0.0010367567652329996, -0.00459519989492062, -0.016877284058530576, 
//				-0.007407790842487786, 5.79801444082945E-4, 0.0010312001807100807, 0.00443890779958875, -0.01694125523540437, 
//				-0.0073708455901976225, 0.0021525019254398483, 0.0010431945012894682, -0.004253044450184699, -0.017074456971173102, 
//				0.00753810565760444, 5.569119120049587E-4, 0.0010453664330827006, 0.004528222381208168, -0.016891046387864164
//		};
		
		//rho = 0.02340
//		double[] parameters = new double[]{
//				0.007501204266326349, 0.0034000358892937767, 4.8766671143592893E-4, -0.003627674501502527, -0.013931054287773574, 
//				-0.006501485514858695, 0.0022123033967697057, 4.81498231578084E-4, 0.0034739426074427556, -0.013991771463094138, 
//				-0.006463658239866977, 0.0034079791057508814, 4.935270294740875E-4, -0.0032873867491024587, -0.014125336723759974, 
//				0.007478521304168808, 0.0021893654117202824, 4.963243715468966E-4, 0.003558280011654309, -0.01394941048125965
//		};

		//rho = 0.02360
//		double[] parameters = new double[]{
//				0.00587858625614348, 3.498136370806539E-4, 2.530802124373233E-4, -0.0024176577048584135, -0.01112289099622431, 
//				-0.006839087612318747, -5.369085724467176E-4, 2.4714008544199687E-4, 0.0022626840421924323, -0.011181846342811854, 
//				-0.0068040941664580775, 3.5690429130867413E-4, 2.6041157724316527E-4, -0.0020879869697364296, -0.011306573172402495, 
//				0.005860486080591986, -5.58272442724777E-4, 2.62622984231827E-4, 0.0023518729234672494, -0.011125642111254046
//		};
		
		//rho = 0.02380
//		double[] parameters = new double[]{
//				0.006280754154976942, 0.0018647826805411045, 2.20680440587407E-4, -0.001400671718729978, -0.00845654317120395, 
//				-0.004716070064603048, 0.0012650164456099395, 2.1514476269008585E-4, 0.001247200935554201, -0.00852300319734449, 
//				-0.00467805245058055, 0.0018722195220994813, 2.2796379324396685E-4, -0.0010672455145830477, -0.00864585101097989, 
//				0.006259429192793395, 0.0012429906588980787, 2.302416130341177E-4, 0.0013344577630835953, -0.008469468776208728
//		};
		
		//rho = 0.02400
//		double[] parameters = new double[]{
//				0.004145316356587442, -9.97188385687702E-4, 0.0010652863572599772, -6.111757291552174E-4, -0.006099678343603974, 
//				-0.004860191100274785, -0.0011184372646839426, 0.0010591048051743549, 4.556229648632361E-4, -0.006156713775670309, 
//				-0.004822918569318898, -9.896000906163874E-4, 0.001072405237441662, -2.7171155047637206E-4, -0.006288187978589316, 
//				0.004124528285182966, -0.0011414286354795309, 0.001074732339805428, 5.412368161856786E-4, -0.006105718687480315
		
//		};
		
		double[] valMin = new double[parameters.length];
		double[] valMax = new double[parameters.length];
	
		double minScale = scale;
		double maxScale = scale;
		
		System.out.println("****************  " + scale + "  *******************");
		boolean reScale = false;
		double reScaleValue = 0.2;
		
		for (int i=0; i<valMin.length; i++){
			if(reScale){
				if(i%20==0||i%20==1||i%20==2||i%20==5||i%20==6||i%20==7
						||i%20==10||i%20==11||i%20==12||i%20==15||i%20==16||i%20==17){
					valMin[i] = parameters[i] - reScaleValue;
				} else {
					valMin[i] = parameters[i] * (1 - minScale);
				}
				
			} else {
				valMin[i] = parameters[i] * (1 - minScale);
			}
		}
		
		for (int i=0; i<valMax.length; i++){
			if(reScale){
				if(i%20==0||i%20==1||i%20==2||i%20==5||i%20==6||i%20==7
						||i%20==10||i%20==11||i%20==12||i%20==15||i%20==16||i%20==17){
					valMax[i] = parameters[i] + reScaleValue;
				} else {
					valMax[i] = parameters[i] * (1 + maxScale);
				}
				
			} else {
				valMax[i] = parameters[i] * (1 + maxScale);
			}
		}
		
//		double[] parameters = new double[]{1.631, 
//				0.00, 0.00, 0.0, 0.0, 0.0, 
//				0.00, 0.00, 0.0, 0.0, 0.0,  
//				0.00, 0.00, 0.0, 0.0, 0.0, 
//				0.00, 0.00, 0.0, 0.0, 0.0};
		
		
		MinimizeBetaNitrogenLatticeParameterFromFile func = new MinimizeBetaNitrogenLatticeParameterFromFile(Space.getInstance(3), nC, density, parameters);
		System.out.println("initial energy: " + func.getEnergy(parameters));

		func.doFindMinimum(valMin, valMax, parameters);

		try {
			FileWriter fileWriter = new FileWriter(filename,false);
			
			for (int i=0; i<parameters.length; i++){
			
				fileWriter.write(parameters[i]+" ");
				
				if(i>0&&(i+1)%5==0){
					fileWriter.write("\n");
						
				}
				//fileWriter.write(parameters[i]+ "\n");
			}
			fileWriter.close();
			
		} catch(IOException e){
			throw new RuntimeException("Failed to write file!!" + e);
		
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
	protected double[] energy;
	protected double[] allValue;
	
	protected double initLat, latticeEnergy;
	protected int numMolecule;
	private static final long serialVersionUID = 1L;
}
