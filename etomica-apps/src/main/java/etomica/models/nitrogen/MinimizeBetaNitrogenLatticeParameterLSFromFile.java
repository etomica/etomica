/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import java.io.FileWriter;
import java.io.IOException;

import etomica.box.Box;
import etomica.atom.IMoleculeList;
import etomica.data.DataInfo;
import etomica.data.FunctionData;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveTriclinic;
import etomica.models.nitrogen.LatticeSumCrystalMolecular.DataGroupLSC;
import etomica.normalmode.BasisBigCell;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.units.Energy;

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
 * with lattice sum
 * 
 * 
 * @author Tai Boon Tan
 *
 */
public class MinimizeBetaNitrogenLatticeParameterLSFromFile extends Simulation {
	
	public MinimizeBetaNitrogenLatticeParameterLSFromFile(Space space, double density, double[] u, double rC){
		super(space);
		this.space = space;
		this.density = density;
		
		double ratio = 1.631;
		double aDim = Math.pow(4.0/(Math.sqrt(3.0)*ratio*density), 1.0/3.0);
		double cDim = aDim*ratio;
		//System.out.println("aDim: " + aDim + " ;cDim: " + cDim);
		
		int [] nCells = new int[]{1,2,1};
		Basis basisHCP = new BasisHcp();
		basis = new BasisBigCell(space, basisHCP, nCells);
        
		ConformationNitrogen conformation = new ConformationNitrogen(space);
		SpeciesN2 species = new SpeciesN2(space);
		species.setConformation(conformation);
		addSpecies(species);
		
		SpeciesN2B ghostSpecies = new SpeciesN2B(space);
		ghostSpecies.setConformation(conformation);
		addSpecies(ghostSpecies);
		
		int numMolecule = 4;
		box = new Box(space);
		addBox(box);
		box.setNMolecules(species, numMolecule);		
		
		ghostBox = new Box(space);
		addBox(ghostBox);
		ghostBox.setNMolecules(ghostSpecies, 1);
		
		primitive = new PrimitiveTriclinic(space, aDim, 2*aDim, cDim, Math.PI*(90/180.0),Math.PI*(90/180.0),Math.PI*(120/180.0));

		double param[][] = new double[4][5];
		for(int i=0; i<4; i++){
			for(int j=0; j<5; j++){
				param[i][j] = u[i*5+j];
					
			}	
		}
		
		coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordinateDef.setIsBetaLatticeSum();
		coordinateDef.setIsDoLatticeSum();
		coordinateDef.setOrientationVectorBetaLatticeSum(space, density, param);
		coordinateDef.initializeCoordinates(new int[]{1,1,1});
		
		potential = new P2Nitrogen(space, rC);
		potential.setBox(box);
		potential.setEnablePBC(false);
		
		this.nLayer = (int)(rC/aDim+0.5);
		
		FunctionData<Object> function = new FunctionData<Object>() {
			public IData f(Object obj) {
				data.x = potential.energy((IMoleculeList)obj);
				return data;
			}
			public IDataInfo getDataInfo() {
				return dataInfo;
			}
			final DataInfo dataInfo = new DataDouble.DataInfoDouble("Lattice energy", Energy.DIMENSION);
			final DataDouble data = new DataDouble();
		};
		
		BravaisLatticeCrystal lattice = new BravaisLatticeCrystal(primitive, basis);
		LatticeSumCrystalMolecular latticeSum = new LatticeSumCrystalMolecular(lattice, coordinateDef, ghostBox);
		latticeSum.setMaxLatticeShell(nLayer);
		
		double sum = 0;
	    double basisDim = lattice.getBasis().getScaledCoordinates().length;
		DataGroupLSC data = (DataGroupLSC)latticeSum.calculateSum(function);
        for(int j=0; j<basisDim; j++) {
            for(int jp=0; jp<basisDim; jp++) {
                sum += ((DataDouble)data.getDataReal(j,jp)).x; 
            }
        }

		System.out.println("initial energy: " + (0.5*sum/basisDim));
	}
	

	public double getEnergy (double[] u){

		double param[][] = new double[4][5];
		for(int i=0; i<4; i++){
			for(int j=0; j<5; j++){
				param[i][j] = u[i*5+j];
					
			}	
		}
		
		CoordinateDefinitionNitrogen coordDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordDef.setIsBetaLatticeSum();
		coordDef.setIsDoLatticeSum();
		coordDef.setOrientationVectorBetaLatticeSum(space, density, param);
		coordDef.initializeCoordinates(new int[]{1,1,1});
		
		FunctionData<Object> function = new FunctionData<Object>() {
			public IData f(Object obj) {
				data.x = potential.energy((IMoleculeList)obj);
				return data;
			}
			public IDataInfo getDataInfo() {
				return dataInfo;
			}
			final DataInfo dataInfo = new DataDouble.DataInfoDouble("Lattice energy", Energy.DIMENSION);
			final DataDouble data = new DataDouble();
		};
		
		BravaisLatticeCrystal lattice = new BravaisLatticeCrystal(primitive, basis);
		LatticeSumCrystalMolecular latticeSum = new LatticeSumCrystalMolecular(lattice, coordDef, ghostBox);
		latticeSum.setMaxLatticeShell(nLayer);
		
		double sum = 0;
	    double basisDim = lattice.getBasis().getScaledCoordinates().length;
		DataGroupLSC data = (DataGroupLSC)latticeSum.calculateSum(function);
        for(int j=0; j<basisDim; j++) {
            for(int jp=0; jp<basisDim; jp++) {
                sum += ((DataDouble)data.getDataReal(j,jp)).x; 
            }
        }

		return 0.5*sum/basisDim;
	}
	

	public void doFindMinimum(double[] minVal, double[] maxVal, double[] parameter){
		this.parameters = parameter;
		int numIter = 1;
		double initEnergy = 0.0;
		double afterEnergy = 0.0;
		double initParam = 0.0;
		
		while(numIter < 50){
			for (int iVar=0; iVar<parameter.length; iVar++){

				initEnergy = getEnergy(parameters);
				initParam = parameters[iVar];
				
				parameters[iVar] = findOptParameter(minVal[iVar], maxVal[iVar], parameters, iVar);
				
				afterEnergy = getEnergy(parameters);
				
				if(afterEnergy < initEnergy){
					System.out.println("**************** LOWER ENERGY LATTICE STRUCTURE FOUND! *******************");
					for(int i=0; i<parameter.length; i++){
						System.out.print(parameters[i]+", ");
						if((i+1)%5==0){
							System.out.println("");
						}
					}
					System.out.println(numIter + " "+ iVar+" lower lattice energy (sim unit): " + getEnergy(parameters));
					
					
//					if(Math.abs(parameters[iVar]) < 1e-8){
//		            	minVal[iVar] = parameter[iVar] - 0.1;
//						maxVal[iVar] = parameter[iVar] + 0.1;
//			        	
//		            } else{
//						minVal[iVar] = parameter[iVar] *(1-0.5/Math.sqrt(numIter));
//						maxVal[iVar] = parameter[iVar] *(1+0.5/Math.sqrt(numIter));
//			        
//		            }
										
				} else {
					parameters[iVar] = initParam;
					
				}
			}
			
			++numIter;
		}
		
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

                if ((energy[1] > energy[0] && energy[1] > energy[2])){ 
                    // we found a maximum, due to numerical precision failure
                    // just bail and pretend that the middle point is the global minimum
                	//System.out.println("first one ");
                    return allValue[1];
                }
                
                if((energy[1] > energy[0] && energy[1] > energy[2])  ){
                	//System.out.println("second one ");
                	return allValue[0];
                }
                
                if((energy[1] < energy[0] && energy[1] > energy[2])  ){
                	//System.out.println("third one ");
                	return allValue[2];
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
                	//System.out.println("***"+value + " ;0: " + allValue[0]+ " ;1: " + allValue[1]+ " ;2: " + allValue[2] );
                	//System.out.println("middle");
                    return value;
                }
                
                if (Math.abs(energy[0]-energy[1])<1e-12 || Math.abs(energy[1]-energy[2])<1e-12 ||Math.abs(energy[0]-energy[2])<1e-12) {
                	
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
                	//System.out.println("last one");
                	return value;
                }
            }
        }
    }
	
	public double getLatticeEnergy(){
		return latticeEnergy;
	}
	
	public static void main(String[] args){
		
		String filename = "/tmp/inputd0.02400";
		double density = 0.0230;
		double scale = 0.001;
		boolean reScale = true;
		double reScaleValue = 0.0058;
		double rC = 100;
		
	    if(args.length > 0){
			filename = args[0];
		}
		if(args.length > 1){
			density = Double.parseDouble(args[1]);
		}
		if(args.length > 2){
			scale = Double.parseDouble(args[3]);
		}
        
//		double[] parameters = new double[21];
//		double[][] paramFromFile = ArrayReader1D.getFromFile(filename);
//		for (int i=0; i<parameters.length;i++){
//			parameters[i] = paramFromFile[i][0];
//		}
		
//		double[] parameters = new double[20];
//		double[][] paramFromFile = ArrayReader1D.getFromFile(filename);
//		int k=0;
//		for (int i=0; i<paramFromFile.length;i++){
//			for (int j=0; j<paramFromFile[0].length;j++){
//				
//				parameters[k]=paramFromFile[i][j];
//				k++;
//			}	
//		}

//		for(int i=0; i<parameters.length; i++){
//			System.out.print(parameters[i]+", ");
//			if((i+1)%5==0){
//				System.out.println("");
//			}
//		}
//		
//		System.exit(1);
		//rho = 0.02300
		double[] parameters = new double[]{		
				//initial energy: -839945.2665296119
//				0.008643467781612603, 0.0026271868875943002+0.02, 0.0019564392850298443, -0.005368542841198517, -0.01963128167700018, 
//				-0.007220434499475046, 7.968809489413192E-4, 0.0019507333095460664, 0.0052101555943571045, -0.019695510442517328, 
//				-0.00718296472671309, 0.002634905247602907, 0.001963151739179857, -0.005024039971913238, -0.019828787213698948, 
//				0.008622065304767177, 7.740172242212951E-4, 0.0019653547372104462, 0.005300512113324435, -0.019643975672095665
				
//				0.008638157179678502, 0.0026096108999515712, 0.001965450797937235, -0.005199377307519452, -0.01949139253429572, 
//				-0.00722492380384562, 7.955824100638133E-4, 0.001959658798801882, 0.005040994127700166, -0.019555913209941954, 
//				-0.007187648130651289, 0.002616949038447896, 0.001971577869489544, -0.0048573388140697666, -0.019690309240519566, 
//				0.008617185563945931, 7.729517745915459E-4, 0.001973779393155908, 0.005133992160722994, -0.01950522916318097, 
//				5 10 lower lattice energy (sim unit): -828.8472214747077 rC = 50A
				
				0.008637993709530622, 0.0026096846786469183, 0.001965450797937235, -0.005199377307519452, -0.019494049784189182, 
				-0.007224721886565305, 7.955824100638133E-4, 0.001959658798801882, 0.0050413936876156145, -0.019558199228975583, 
				-0.007187766723770279, 0.002617189177234774, 0.0019716996369754576, -0.004856835391208816, -0.019691546165277733, 
				0.008617256066040527, 7.728350550139484E-4, 0.00197391904018354, 0.005133136709939338, -0.019507046532061444, 
//				1 17 lower lattice energy (sim unit): -829.0483133535762
		};
		
		//rho = 0.02320
//		double[] parameters = new double[]{
//				//initial energy: -838468.2470407118
//				0.007559425210441226, 0.002144669727082795, 0.0010366834714885741, -0.004595884623994312, -0.016877284058530576, 
//				-0.007407790842487786, 5.798327761478291E-4, 0.0010309478744233874, 0.004438994439937659, -0.01694125523540437, 
//				-0.0073708455901976225, 0.0021522528745652747, 0.0010431945012894682, -0.0042530444501846995, -0.017074456971173102, 
//				0.007538105657604441, 5.570274852839988E-4, 0.0010453664330827006, 0.004528222381208168, -0.016891046387864164
//		};
		
		//rho = 0.02340
//		double[] parameters = new double[]{
//				//initial energy: -836355.0136543293
//				0.0075027058787762516, 0.0033989513550305036, 4.8687083861097063E-4, -0.0036299483934711984, -0.01392707117422452, 
//				-0.006498710169229157, 0.00221226100907468, 4.811043882164454E-4, 0.003471764164319288, -0.013991813503532121, 
//				-0.0064612608166665095, 0.0034066524637804117, 4.934453566205766E-4, -0.0032860156263545826, -0.014124884852503673, 
//				0.007481336660995784, 0.002189447294110959, 4.956782260127964E-4, 0.0035610933667526786, -0.013940641137692313
//				
//		};

		//rho = 0.02360
//		double[] parameters = new double[]{
//				//initial energy: -833597.8814881488
//				0.00587985895541576, 3.318719115741678E-4, 2.590043720731268E-4, -0.0024202818210948412, -0.011108812803056978, 
//				-0.006840799353167392, -5.522148544040112E-4, 2.5340620457769037E-4, 0.0022629703361688575, -0.01117312661414254, 
//				-0.0068037096583788426, 3.3965397522564555E-4, 2.6566579337048335E-4, -0.002077526073115998, -0.011306619503413384, 
//				0.0058585116614836355, -5.7505802693964E-4, 2.677759602717458E-4, 0.0023510309480724322, -0.011122477252704037
//				
//		};
		
		//rho = 0.02380
//		double[] parameters = new double[]{
//				//initial energy: -830194.1281584246
//				0.006281296090595155, 0.0018636960324207065, 2.2110240804828693E-4, -0.0014036809054544244, -0.008449971754047568, 
//				-0.004714745613562198, 0.0012650641596351473, 2.1536173401871355E-4, 0.0012467939175371222, -0.008514524358720084, 
//				-0.004677708017512473, 0.001871373475529715, 2.275760062626014E-4, -0.0010607862675041499, -0.00864766987910912, 
//				0.006259751637666007, 0.0012419882772394947, 2.299559332681835E-4, 0.0013340679895082646, -0.008464163570206894
//		};
		
//		//rho = 0.02400
//		double[] parameters = new double[]{
//				//initial energy: -826145.2631488933
//				0.0041455759544500636, -9.973125892194278E-4, 0.0010654667332234066, -6.112646497823689E-4, -0.006091543797628076, 
//				-0.004860191100274781, -0.0011185047819778846, 0.0010598981530198438, 4.553686694383018E-4, -0.006154751297776442, 
//				-0.004822876003367027, -9.896044425178226E-4, 0.0010722741750401753, -2.699940955128142E-4, -0.006287800037644388, 
//				0.00412436795091807, -0.0011411936940946507, 0.001074384513588608, 5.421605014910879E-4, -0.0061044675703691605
//		};
		
		double[] valMin = new double[parameters.length];
		double[] valMax = new double[parameters.length];
	
		double minScale = scale*1.01;
		double maxScale = scale*1.02;
		
		System.out.println("****************  " + scale + "  *******************");

		for (int i=0; i<valMin.length; i++){
			if(reScale){
				//if(i%20==0||i%20==1||i%20==2||i%20==5||i%20==6||i%20==7
				//		||i%20==10||i%20==11||i%20==12||i%20==15||i%20==16||i%20==17){
					valMin[i] = parameters[i] - reScaleValue*0.99;
				//} 
//				else {
//					valMin[i] = parameters[i] * (1 - minScale);
//				}
				
			} else {
				valMin[i] = parameters[i] * (1 - minScale);
			}
		}
		
		for (int i=0; i<valMax.length; i++){
			if(reScale){
				//if(i%20==0||i%20==1||i%20==2||i%20==5||i%20==6||i%20==7
				//		||i%20==10||i%20==11||i%20==12||i%20==15||i%20==16||i%20==17){
					valMax[i] = parameters[i] + reScaleValue;
				//} 
//				else {
//					valMax[i] = parameters[i] * (1 + maxScale);
//				}
				
			} else {
				valMax[i] = parameters[i] * (1 + maxScale);
			}
		}
		
		
		MinimizeBetaNitrogenLatticeParameterLSFromFile func = new MinimizeBetaNitrogenLatticeParameterLSFromFile(Space.getInstance(3), density, parameters, rC);
		
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
	
	protected CoordinateDefinitionNitrogen coordinateDef;
	protected P2Nitrogen potential;
	protected Box box, ghostBox;
	protected SpeciesN2 species;
	protected double density;
	protected Basis basis;
	protected Space space;
	protected Primitive primitive;
	protected double [] parameters;  
	protected int[] nC;
	protected double[] energy;
	protected double[] allValue;
	protected int nLayer;
	
	protected double latticeEnergy;
	private static final long serialVersionUID = 1L;
}
