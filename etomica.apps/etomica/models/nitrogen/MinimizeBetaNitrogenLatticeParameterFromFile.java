package etomica.models.nitrogen;

import java.io.FileWriter;
import java.io.IOException;

import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.normalmode.BasisBigCell;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.units.Degree;
import etomica.util.numerical.ArrayReader1D;

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
	

	public MinimizeBetaNitrogenLatticeParameterFromFile(ISpace space, int[] nC, double density, double[] u){
		super(space);
		this.space = space;
		this.density = density;
		this.nC = nC;
		
		double ratio = 1.631;
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
		
		boxDim = new IVector[3];
		boxDim[0] = space.makeVector(new double[]{nC[0]*aDim, 0, 0});
		boxDim[1] = space.makeVector(new double[]{-nC[1]*aDim*Math.cos(Degree.UNIT.toSim(60)), nC[1]*aDim*Math.sin(Degree.UNIT.toSim(60)), 0});
		boxDim[2] = space.makeVector(new double[]{0, 0, nC[2]*cDim});
		
		boundary = new BoundaryDeformablePeriodic(space, boxDim);
		primitive = new PrimitiveHexagonal(space, nC[0]*aDim, nC[2]*cDim);
		
		coordinateDefinition = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordinateDefinition.setIsBeta();
		coordinateDefinition.setOrientationVectorBeta(space);
		coordinateDefinition.initializeCoordinates(new int[]{1,1,1});
		
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
		
		box.setBoundary(boundary);
		double rC = nC[0]*aDim*0.475;
		//System.out.println("rC: " + rC);
		P2Nitrogen potential = new P2Nitrogen(space, rC);
		potential.setBox(box);
		
		potentialMaster.addPotential(potential, new ISpecies[]{species, species});
		
		meterPotential = new MeterPotentialEnergy(potentialMaster);
    	meterPotential.setBox(box);
    	initLat = meterPotential.getDataAsScalar()/numMolecule;
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
		
		while(numIter < 30){
		
			for (int iVar=1; iVar<parameter.length; iVar++){
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
					//System.out.println("**************** LOWER ENERGY LATTICE STRUCTURE FOUND! *******************");
					System.out.println(numIter + " "+ iVar+" lattice energy (sim unit): " + getEnergy(parameters)/numMolecule);
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
		
		String filename = "/tmp/inputd0.021";
		double density = 0.021;
		int nCells = 8;
		double scale = 1.2;
		
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
		double[] parameters = new double[21];
		double[][] paramFromFile = ArrayReader1D.getFromFile(filename);
		for (int i=0; i<parameters.length;i++){
			parameters[i] = paramFromFile[i][0];
		}
		
//		double[] parameters = new double[]{1.631, 
//				0.00967627511673874, -4.0335087762077794E-4, -4.979390112773919E-4, -0.00621165210733723, -0.026684242064978696,
//				 -0.010964121568701888, -8.89076558191218E-4, -5.103540297036002E-4, 0.006075743547189045, -0.02662786102380948,
//				 -0.010897042866236585, -3.862882969924215E-4, -4.957174274794356E-4, -0.00589031685970485, -0.026762492129711716,
//				 0.00965324315670293, -9.269032792479676E-4, -4.8760886758373234E-4, 0.006107404986730156, -0.026621244306518587
//
//		};

		
		
		double[] valMin = new double[parameters.length];
		double[] valMax = new double[parameters.length];
		
		double minScale = scale;
		double maxScale = scale;
		
		System.out.println("****************  " + scale + "  *******************");
		boolean reScale = false;
		double reScaleValue = 0.06;
		for (int i=1; i<valMin.length; i++){
			if(reScale){
				if(i%20==1||i%20==2||i%20==3||i%20==6||i%20==7||i%20==8
						||i%20==11||i%20==12||i%20==13||i%20==16||i%20==17||i%20==18){
					valMin[i] = parameters[i] - reScaleValue;
				} else {
					valMin[i] = parameters[i] * (1 - minScale);
				}
				
			} else {
				valMin[i] = parameters[i] * (1 - minScale);
			}
		}
		
		for (int i=1; i<valMax.length; i++){
			if(reScale){
				if(i%20==1||i%20==2||i%20==3||i%20==6||i%20==7||i%20==8
						||i%20==11||i%20==12||i%20==13||i%20==16||i%20==17||i%20==18){
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
		func.doFindMinimum(valMin, valMax, parameters);

		try {
			FileWriter fileWriter = new FileWriter(filename,false);
			
			for (int i=0; i<parameters.length; i++){
			
				fileWriter.write(parameters[i]+ "\n");
			}
			fileWriter.write(func.getEnergy(parameters)/func.numMolecule + "\n");
			fileWriter.close();
			
		} catch(IOException e){
			throw new RuntimeException("Failed to write coord data normalize coord U" + e);
		
		}
	}
	
	protected CoordinateDefinitionNitrogen coordinateDefinition;
	protected MeterPotentialEnergy meterPotential;
	protected PotentialMaster potentialMaster;
	protected IBox box;
	protected SpeciesN2 species;
	protected double density;
	protected ISpace space;
	protected Primitive primitive;
	protected Boundary boundary;
	protected IVector[] boxDim; 
	protected double aDim, cDim;
	protected double [] parameters;  
	protected int[] nC;
	protected double[] energy;
	protected double[] allValue;
	
	protected double initLat, latticeEnergy;
	protected int numMolecule;
	private static final long serialVersionUID = 1L;
}
