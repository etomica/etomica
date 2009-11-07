package etomica.normalmode;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorAverageFixedOutputFile;
import etomica.data.DataPump;
import etomica.data.IEtomicaDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.listener.IntegratorListenerAction;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;

/**
 * Simulation to sample Umbrella Overlap Sampling
 * 
 * The original Umbrella Sampling Simulation
 * 	- used to check for the computation time
 * 
 * @author Tai Boon Tan
 */
public class SimUmbrellaSoftSphere extends Simulation {

	private static final String APP_NAME = "Sim Umbrella's";

    public SimUmbrellaSoftSphere(Space _space, int numAtoms, double density, double temperature, String filename, int exponent) {
        super(_space);

        String refFileName = filename +"_ref";
        FileReader refFileReader;
        try {
        	refFileReader = new FileReader(refFileName);
        } catch (IOException e){
        	throw new RuntimeException ("Cannot find refPref file!! "+e.getMessage() );
        }
        
        try {
        	BufferedReader bufReader = new BufferedReader(refFileReader);
        	String line = bufReader.readLine();
        	
        	refPref = Double.parseDouble(line);
        	setRefPref(refPref);
        	
        } catch (IOException e){
        	throw new RuntimeException(" Cannot read from file "+ refFileName);
        }
        
        int D = space.D();
        
        BasisNMOneCell scaledBasis = new BasisNMOneCell(_space, numAtoms, density);
                
        potentialMasterMonatomic = new PotentialMasterMonatomic(this);
        integrator = new IntegratorMC(potentialMasterMonatomic, getRandom(), temperature);
       
        species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        //Target        
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);

        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
      
       	double L = Math.pow(4.0/density, 1.0/3.0);
       	int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
        primitive = new PrimitiveCubic(space, n*L);
        nCells = new int[]{1,1,1};
        boundary = new BoundaryRectangularPeriodic(space, n*L);
        basis = new Basis(scaledBasis.getScaledBasis());
        
        box.setBoundary(boundary);
        
        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        normalModes = new NormalModesFromFile(filename, D);
        /*
         * nuke this line when it is derivative-based
         */
        normalModes.setTemperature(temperature);
        coordinateDefinition.initializeCoordinates(nCells);
        
        Potential2SoftSpherical potential = new P2SoftSphere(space, 1.0, 1.0, exponent);
        double truncationRadius = boundary.getBoxSize().getX(0) * 0.495;
        P2SoftSphericalTruncatedShifted pTruncated = new P2SoftSphericalTruncatedShifted(space, potential, truncationRadius);
        IAtomType sphereType = species.getLeafType();
        potentialMasterMonatomic.addPotential(pTruncated, new IAtomType[] { sphereType, sphereType });
        
        integrator.setBox(box);
        
        /*
         *  1-body Potential to Constraint the atom from moving too far
         *  	away from its lattice-site
         *  
         */

        P1Constraint p1Constraint = new P1Constraint(space, primitive, box, coordinateDefinition);
        potentialMasterMonatomic.addPotential(p1Constraint, new IAtomType[] {sphereType});
        
        potentialMasterMonatomic.lrcMaster().setEnabled(false);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMasterMonatomic);
        meterPE.setBox(box);
        latticeEnergy = meterPE.getDataAsScalar();
        
        move = new MCMoveAtomCoupledUmbrella(potentialMasterMonatomic, getRandom(), 
        		coordinateDefinition, normalModes, getRefPref(), space);
        move.setTemperature(temperature);
        move.setLatticeEnergy(latticeEnergy);
        integrator.getMoveManager().addMCMove(move);
        ((MCMoveStepTracker)move.getTracker()).setNoisyAdjustment(true);
      
        meterHarmonicEnergy = new MeterHarmonicEnergy(coordinateDefinition, normalModes);
        
        meterEnergy = new MeterPotentialEnergy(potentialMasterMonatomic);
        meterEnergy.setBox(box);
        
    }


	public double getRefPref() {
		return refPref;
	}

	public void setRefPref(double refPref) {
		this.refPref = refPref;
	}
    
    /**
     * @param args
     */
    public static void main(String[] args) {
      	
        SimBennetParam params = new SimBennetParam();
        String inputFilename = null;
        if (args.length > 0) {
            inputFilename = args[0];
        }
        if (inputFilename != null) {
            ReadParameters readParameters = new ReadParameters(inputFilename, params);
            readParameters.readParameters();
        }
        double density = params.density/1000;
        int exponentN = params.exponentN;
        long numSteps = params.numSteps;
        final int numAtoms = params.numMolecules;
        double temperature = params.temperature;
        double harmonicFudge = params.harmonicFudge;
        int D = params.D;
        String filename = params.filename;
        if (filename.length() == 0) {
        	System.err.println("Need input files!!!");
            filename = "FCC_SoftSphere_n"+exponentN+"_T"+ (int)Math.round(temperature*10);
        }
        
        final long startTime = System.currentTimeMillis();
        System.out.println("Start Time: " + startTime);
        System.out.println("Running "+(D==1 ? "1D" : (D==3 ? "FCC" : "2D hexagonal")) +" soft sphere Umbrella's Sampling simulation");
        System.out.println(numAtoms+" atoms at density "+density+" and temperature "+temperature);
        System.out.println("exponent N: "+ exponentN );
        System.out.println("total steps: "+ numSteps);
        System.out.println("output data to "+filename);
        
        //construct simulation
        final SimUmbrellaSoftSphere sim = new SimUmbrellaSoftSphere(Space.getInstance(D), numAtoms, density, temperature, filename, exponentN);
        
        IEtomicaDataSource[] samplingMeters = new IEtomicaDataSource[2];
        
        sim.activityIntegrate.setMaxSteps(numSteps/10);
        sim.getController().actionPerformed();
        System.out.println("System Equilibrated!");
        
        sim.getController().reset();

        /*
         * 
         */
        
        // Harmonic Sampling
        
        final MeterSamplingHarmonic meterSamplingHarmonic = new MeterSamplingHarmonic(sim.integrator, sim.move);
        meterSamplingHarmonic.setRefPref(sim.refPref);
        samplingMeters[0] = meterSamplingHarmonic;
        
//        final AccumulatorAverageFixedOutputFile dataAverageSamplingHarmonic = new AccumulatorAverageFixedOutputFile();
//        dataAverageSamplingHarmonic.setFile(filename+"_UmblnQharm");
        
        final AccumulatorAverageFixed dataAverageSamplingHarmonic = new AccumulatorAverageFixed();
        
        DataPump pumpSamplingHarmonic = new DataPump(samplingMeters[0], dataAverageSamplingHarmonic);
        dataAverageSamplingHarmonic.setBlockSize(200);
        IntegratorListenerAction pumpSamplingHarmonicListener = new IntegratorListenerAction(pumpSamplingHarmonic);
        sim.integrator.getEventManager().addListener(pumpSamplingHarmonicListener);
    
        
        // Target Sampling
        
        final MeterSamplingTarget meterSamplingTarget = new MeterSamplingTarget(sim.integrator, sim.move);
        meterSamplingTarget.setRefPref(sim.refPref);
        samplingMeters[1] = meterSamplingTarget;
        
//        final AccumulatorAverageFixedOutputFile dataAverageSamplingTarget = new AccumulatorAverageFixedOutputFile();
//        dataAverageSamplingTarget.setFile(filename+"_UmblnQtarg");
        
        final AccumulatorAverageFixed dataAverageSamplingTarget = new AccumulatorAverageFixed();
        
        DataPump pumpSamplingTarget = new DataPump(samplingMeters[1], dataAverageSamplingTarget);
        dataAverageSamplingTarget.setBlockSize(200);
        IntegratorListenerAction pumpSamplingTargetListener = new IntegratorListenerAction(pumpSamplingTarget);
        sim.integrator.getEventManager().addListener(pumpSamplingTargetListener);
        
        if (numAtoms == 32){
            pumpSamplingHarmonicListener.setInterval(100);
            pumpSamplingTargetListener.setInterval(100);
        } else if (numAtoms >= 108){
        	 pumpSamplingHarmonicListener.setInterval(300);
            pumpSamplingTargetListener.setInterval(300);
            /*
            if (temperature >= 1.1){
            	dataAverageSamplingHarmonic.setBlockSize(200);
            	dataAverageSamplingTarget.setBlockSize(200);
            }*/
        } 
        
        /*else {
            pumpSamplingHarmonicListener.setInterval(1);
            pumpSamplingTargetListener.setInterval(1);
        }
       */
        
		double[][] omega2 = sim.normalModes.getOmegaSquared();
		double[] coeffs = sim.normalModes.getWaveVectorFactory().getCoefficients();
		double AHarmonic = 0;
		for(int i=0; i<omega2.length; i++) {
		      for(int j=0; j<omega2[0].length; j++) {
		            if (!Double.isInfinite(omega2[i][j])) {
		               AHarmonic += coeffs[i]*Math.log(omega2[i][j]*coeffs[i]/(temperature*Math.PI));
		            }
		      }
		}

        int totalCells = 1;
        for (int i=0; i<D; i++) {
        		totalCells *= sim.nCells[i];
        }
		int basisSize = sim.basis.getScaledCoordinates().length;
		double fac = 1;
		if (totalCells % 2 == 0) {
			fac = Math.pow(2,D);
		}
		AHarmonic -= Math.log(Math.pow(2.0, basisSize*D*(totalCells - fac)/2.0)); // Math.pow(totalCells,0.5*D));
		System.out.println("Harmonic-reference free energy: "+AHarmonic*temperature + " " +AHarmonic*temperature/(numAtoms));
		System.out.println(" ");
		
//		FileWriter fileWriterHarm, fileWriterTarget;
//		
//		try{
//			fileWriterHarm = new FileWriter(filename+"_Qharm");
//			fileWriterTarget = new FileWriter(filename+"_Qtarg");
//			
//		} catch (IOException e){
//			fileWriterHarm = null;
//			fileWriterTarget = null;
//		}
//		
//        final double temp = temperature;
//        final double AHarm = AHarmonic;
//        
//        final FileWriter fHarm = fileWriterHarm;
//        final FileWriter fTarg = fileWriterTarget;
//        
//		IAction output = new IAction(){
//			public void actionPerformed(){
//				/*
//				 * Qharmonic = < e0 / [sqrt(e1^2 + alpha^2 * e0^2)]>umbrella
//				 *  Qtarget  = < e1 / [sqrt(e1^2 + alpha^2 * e0^2)]>umbrella
//				 * 
//				 */
//				long currentStep = sim.integrator.getStepCount();
//	
//				double aveQharmonic = dataAverageSamplingHarmonic.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
//				double aveQtarget = dataAverageSamplingTarget.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
//				
//				double Qharmonic = meterSamplingHarmonic.getData().getValue(0);
//				double Qtarget = meterSamplingTarget.getData().getValue(0);
//				
//			    double eQharmonic = dataAverageSamplingHarmonic.getData().getValue(AccumulatorAverage.StatType.ERROR.index);
//			    double eQtarget = dataAverageSamplingTarget.getData().getValue(AccumulatorAverage.StatType.ERROR.index);
//				
//			    /*
//				 * deltaFE_harmonic: beta*(FE_harmonic - FE_umbrella) = - ln(Qharmonic)
//				 *  deltaFE_target : beta*(FE_target - FE_umbrella) = - ln(Qtarget)
//				 */
//				double deltaFE_harmonic = - Math.log(aveQharmonic);
//				double deltaFE_target = - Math.log(aveQtarget);
//				double deltaFE = temp*deltaFE_target - temp*deltaFE_harmonic;
//        
//				long currentTime = System.currentTimeMillis();
//				double error = temp*Math.sqrt( (eQharmonic/aveQharmonic)*(eQharmonic/aveQharmonic) 
//						+ (eQtarget/aveQtarget)*(eQtarget/aveQtarget))/numAtoms;
//				/*
//				if (Double.isNaN(error)){
//					System.out.println("error: "+ error);
//				}
//				System.out.println(currentStep+ " : Time: " + (currentTime - startTime)+ 
//						" ,Targ_FE/N: "+(temp*AHarm+ deltaFE)/numAtoms +
//						" ,error: "+error+
//						" ,Qharmonic: " + aveQharmonic + " , " + eQharmonic +
//						" ,Qtarget: " + aveQtarget + " , "+ eQtarget);
//				*/
//				try {	
//					fHarm.write(currentStep + " " + Qharmonic + " " + aveQharmonic + " " + eQharmonic + "\n");
//					fTarg.write(currentStep + " " + Qtarget + " " + aveQtarget + " " + eQtarget + "\n");
//					
//				} catch (IOException e){
//					
//				}
//				
//			}
//		};
//		
//		IntegratorListenerAction outputListener = new IntegratorListenerAction(output);
//		outputListener.setInterval((int)numSteps/200);
//		sim.integrator.getEventManager().addListener(outputListener);
	    
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.getController().actionPerformed();
        
//		dataAverageSamplingHarmonic.closeFile();
//		dataAverageSamplingTarget.closeFile();
		/*
		 * Qharmonic = < e0 / [e1 + alpha * e0)]>umbrella
		 *  Qtarget  = < e1 / [e1 + alpha * e0)]>umbrella
		 * 
		 */
	    double Qharmonic = dataAverageSamplingHarmonic.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
	    double Qtarget = dataAverageSamplingTarget.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
			
	    double eQharmonic = dataAverageSamplingHarmonic.getData().getValue(AccumulatorAverage.StatType.ERROR.index);
	    double eQtarget = dataAverageSamplingTarget.getData().getValue(AccumulatorAverage.StatType.ERROR.index);
			        
	    /*
	     * deltaFE_harmonic: beta*(FE_harmonic - FE_umbrella) = - ln(Qharmonic)
	     *  deltaFE_target : beta*(FE_target - FE_umbrella) = - ln(Qtarget)
	     */
	    double deltaFE_harmonic = - Math.log(Qharmonic);
	    double deltaFE_target = - Math.log(Qtarget);
	    double deltaFE = temperature*deltaFE_target - temperature*deltaFE_harmonic;
	      
	    System.out.println("Qharmonic: "+ Qharmonic+ " ;error: " +eQharmonic);
	    System.out.println("Qtarget: "+ Qtarget+ " ;error: " +eQtarget);
	    System.out.println("Q: " +(Qtarget/Qharmonic) + " ;error: "+ (Qtarget/Qharmonic)*Math.sqrt((eQharmonic/Qharmonic)*(eQharmonic/Qharmonic) + (eQtarget/Qtarget)*(eQtarget/Qtarget)));
	    System.out.println("\nfree energy difference (umbrella --> harmonic): "+ (temperature*deltaFE_harmonic) + "; error: " + temperature*eQharmonic/Qharmonic);
	    System.out.println("free energy difference (umbrella --> target): "+ (temperature*deltaFE_target) + "; error: " + temperature*eQtarget/Qtarget);
	    	        
	    System.out.println(" beta*deltaFE (harmonic): " + deltaFE_harmonic);
	    System.out.println(" beta*deltaFE (target): " + deltaFE_target);
			        
	    System.out.println("\nHarmonic-reference free energy: " + temperature*AHarmonic);
	    System.out.println("free energy difference (harmonic --> target): " + deltaFE 
	    		           	+ ", error: " + temperature*Math.sqrt( (eQharmonic/Qharmonic)*(eQharmonic/Qharmonic) + (eQtarget/Qtarget)*(eQtarget/Qtarget) ));
		System.out.println("target free energy: " + (temperature*AHarmonic+ deltaFE));
		System.out.println("target free energy per particle: " + (temperature*AHarmonic+ deltaFE)/numAtoms);
		long endTime = System.currentTimeMillis();
		System.out.println("End Time: " + endTime);
		System.out.println("Total time taken: " + (endTime - startTime));
		
//		try {
//			fileWriterHarm.close();
//			fileWriterTarget.close();
//			
//		} catch (IOException e){
//			
//		}
		
    }

    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public IBox box;
    public Boundary boundary;
    public Basis basis;
    public SpeciesSpheresMono species;
    public NormalModes normalModes;
    public int[] nCells;
    public CoordinateDefinition coordinateDefinition;
    public Primitive primitive;
    public PotentialMasterMonatomic potentialMasterMonatomic;
    public double latticeEnergy;
    public MeterHarmonicEnergy meterHarmonicEnergy;
    public MeterPotentialEnergy meterEnergy;
    public MCMoveAtomCoupledUmbrella move;
    public double refPref;
    
    public static class SimBennetParam extends ParameterBase {
    	public int numMolecules = 32;
    	public double density = 1256;
    	public int exponentN = 12;
    	public int D = 3;
    	public long numSteps = 1000000;
    	public double harmonicFudge =1;
    	public String filename = "input_hessian32CBT01";
    	public double temperature = 0.1;
    }

}