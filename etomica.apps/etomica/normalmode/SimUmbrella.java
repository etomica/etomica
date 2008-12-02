package etomica.normalmode;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAction;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistogram;
import etomica.data.DataFork;
import etomica.data.DataLogger;
import etomica.data.DataPump;
import etomica.data.IEtomicaDataSource;
import etomica.data.DataTableWriter;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.DoubleRange;
import etomica.util.HistogramExpanding;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;

/**
 * Simulation to sample Umbrella Sampling
 */
public class SimUmbrella extends Simulation {

	private static final String APP_NAME = "Sim Umbrella's";

    public SimUmbrella(Space _space, int numAtoms, double density, double temperature, String filename, int exponent) {
        super(_space, true);

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
        //System.out.println("refPref is: "+ refPref);
        
        
        int D = space.D();
        
        potentialMasterMonatomic = new PotentialMasterMonatomic(this,space);
        integrator = new IntegratorMC(potentialMasterMonatomic, getRandom(), temperature);
       
        species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);

        //Target        
        box = new Box(this, space);
        addBox(box);
        box.setNMolecules(species, numAtoms);

        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
      
       	double L = Math.pow(4.0/density, 1.0/3.0);
        primitive = new PrimitiveCubic(space, L);
        int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
        nCells = new int[]{n,n,n};
        boundary = new BoundaryRectangularPeriodic(space, getRandom(), n*L);
        basis = new BasisCubicFcc();
        
        box.setBoundary(boundary);
        
        coordinateDefinition = new CoordinateDefinitionLeaf(this, box, primitive, basis, space);
        normalModes = new NormalModesFromFile(filename, D);
        /*
         * nuke this line when it is derivative-based
         */
        //normalModes.setTemperature(temperature);
        coordinateDefinition.initializeCoordinates(nCells);
        
        Potential2SoftSpherical potential = new P2SoftSphere(space);
        double truncationRadius = boundary.getDimensions().x(0) * 0.495;
        P2SoftSphericalTruncatedShifted pTruncated = new P2SoftSphericalTruncatedShifted(space, potential, truncationRadius);
        IAtomTypeLeaf sphereType = species.getLeafType();
        potentialMasterMonatomic.addPotential(pTruncated, new IAtomTypeLeaf[] { sphereType, sphereType });
        
        integrator.setBox(box);
        
        /*
         *  1-body Potential to Constraint the atom from moving too far
         *  	away from its lattice-site
         *  
         */

        P1Constraint p1Constraint = new P1Constraint(space, primitive, box, coordinateDefinition);
        potentialMasterMonatomic.addPotential(p1Constraint, new IAtomTypeLeaf[] {sphereType});
        
        
        potentialMasterMonatomic.lrcMaster().setEnabled(false);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMasterMonatomic);
        meterPE.setBox(box);
        latticeEnergy = meterPE.getDataAsScalar();
        
        MCMoveAtomCoupledUmbrella move = new MCMoveAtomCoupledUmbrella(potentialMasterMonatomic, getRandom(), 
        		coordinateDefinition, normalModes, getRefPref(), space);
        move.setTemperature(temperature);
        move.setLatticeEnergy(latticeEnergy);
        integrator.getMoveManager().addMCMove(move);
      
        meterHarmonicEnergy = new MeterHarmonicEnergy(coordinateDefinition, normalModes);
        meterHarmonicEnergy.setBox(box);
        
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
        int numAtoms = params.numMolecules;
        double temperature = params.temperature;
        double harmonicFudge = params.harmonicFudge;
        int D = params.D;
        String filename = params.filename;
        if (filename.length() == 0) {
        	System.err.println("Need input files!!!");
            filename = "FCC_SoftSphere_n"+exponentN+"_T"+ (int)Math.round(temperature*10);
        }
       
    	
        System.out.println("Running "+(D==1 ? "1D" : (D==3 ? "FCC" : "2D hexagonal")) +" soft sphere Umbrella's Sampling perturbation simulation");
        System.out.println(numAtoms+" atoms at density "+density+" and temperature "+temperature);
        System.out.println("exponent N: "+ exponentN );
        System.out.println("total steps: "+ numSteps);
        System.out.println("output data to "+filename);
        
        //construct simulation
        final SimUmbrella sim = new SimUmbrella(Space.getInstance(D), numAtoms, density, temperature, filename, exponentN);
        
        IEtomicaDataSource[] workMeters = new IEtomicaDataSource[2];
        IEtomicaDataSource[] samplingMeters = new IEtomicaDataSource[2];
        
      // Umbrella Sampling ---> Harmonic
        final MeterWorkUmbrellaHarmonic meterWorkUmbrellaHarmonic = new MeterWorkUmbrellaHarmonic(sim.integrator, sim.meterHarmonicEnergy);
        meterWorkUmbrellaHarmonic.setRefPref(sim.refPref);
        meterWorkUmbrellaHarmonic.setLatticeEnergy(sim.latticeEnergy);
        workMeters[0] = meterWorkUmbrellaHarmonic;
        
        DataFork dataForkHarmonic = new DataFork();
        DataPump pumpHarmonic = new DataPump(workMeters[0], dataForkHarmonic);
        
        final AccumulatorAverageFixed dataAverageHarmonic = new AccumulatorAverageFixed();
        dataForkHarmonic.addDataSink(dataAverageHarmonic);
        sim.integrator.addIntervalAction(pumpHarmonic);
        sim.integrator.setActionInterval(pumpHarmonic, 1);
        
        //Histogram Harmonic
        final AccumulatorHistogram histogramHarmonic = new AccumulatorHistogram(new HistogramExpanding(1, new DoubleRange(0,1)));
        dataForkHarmonic.addDataSink(histogramHarmonic);
        
      // Umbrella Sampling ---> Target
        final MeterWorkUmbrellaTarget meterWorkUmbrellaTarget = new MeterWorkUmbrellaTarget(sim.integrator, sim.meterHarmonicEnergy);
        meterWorkUmbrellaTarget.setRefPref(sim.refPref);
        meterWorkUmbrellaTarget.setLatticeEnergy(sim.latticeEnergy);
        workMeters[1] = meterWorkUmbrellaTarget;
        
        DataFork dataForkTarget = new DataFork();
        DataPump pumpTarget = new DataPump(workMeters[1], dataForkTarget);
        
        final AccumulatorAverageFixed dataAverageTarget = new AccumulatorAverageFixed();
        dataForkTarget.addDataSink(dataAverageTarget);
        sim.integrator.addIntervalAction(pumpTarget);
        sim.integrator.setActionInterval(pumpTarget, numAtoms*2);
        
        //Histogram Target
        final AccumulatorHistogram histogramTarget = new AccumulatorHistogram(new HistogramExpanding(1,new DoubleRange(0,1)));
        dataForkTarget.addDataSink(histogramTarget);

        /*
         * 
         */
        
        // Harmonic Sampling
        
        final MeterSamplingHarmonic meterSamplingHarmonic = new MeterSamplingHarmonic(sim.integrator, sim.meterHarmonicEnergy);
        meterSamplingHarmonic.setRefPref(sim.refPref);
        meterSamplingHarmonic.setLatticeEnergy(sim.latticeEnergy);
        samplingMeters[0] = meterSamplingHarmonic;
        
        final AccumulatorAverageFixed dataAverageSamplingHarmonic = new AccumulatorAverageFixed();
        DataPump pumpSamplingHarmonic = new DataPump(samplingMeters[0], dataAverageSamplingHarmonic);
        sim.integrator.addIntervalAction(pumpSamplingHarmonic);
        sim.integrator.setActionInterval(pumpSamplingHarmonic, 1);
        
        // Target Sampling
        
        final MeterSamplingTarget meterSamplingTarget = new MeterSamplingTarget(sim.integrator, sim.meterHarmonicEnergy);
        meterSamplingTarget.setRefPref(sim.refPref);
        meterSamplingTarget.setLatticeEnergy(sim.latticeEnergy);
        samplingMeters[1] = meterSamplingTarget;
        
        final AccumulatorAverageFixed dataAverageSamplingTarget = new AccumulatorAverageFixed();
        DataPump pumpSamplingTarget = new DataPump(samplingMeters[1], dataAverageSamplingTarget);
        sim.integrator.addIntervalAction(pumpSamplingTarget);
        sim.integrator.setActionInterval(pumpSamplingTarget, 1);
        
        FileWriter fileWriter;
        
        try{
        	fileWriter = new FileWriter(filename+"_SimUmb");
        } catch(IOException e){
        	fileWriter = null;
        }
        
        final String outFileName = filename;
        final double temp = temperature;
        final long steps = numSteps;
        final FileWriter fileWriterUmb = fileWriter;
        
        
        IAction outputAction = new IAction(){
        	public void actionPerformed(){
        		int idStep = sim.integrator.getStepCount();
		        /*
		         * Histogram
		         */
		        //Target
				DataLogger dataLogger = new DataLogger();
				DataTableWriter dataTableWriter = new DataTableWriter();
				dataLogger.setFileName(outFileName + "_hist_UmbTarg");
				dataLogger.setDataSink(dataTableWriter);
				dataTableWriter.setIncludeHeader(false);
				dataLogger.putDataInfo(histogramTarget.getDataInfo());
				
				dataLogger.setWriteInterval(1);
				dataLogger.setAppending(false); //overwrite the file 8/5/08
				dataLogger.putData(histogramTarget.getData());
				dataLogger.closeFile();
		        
				//Harmonic
		        dataLogger.setFileName(outFileName + "_hist_UmbHarm");
		        dataLogger.setAppending(false); // overwrite the file 8/5/08
		        dataLogger.putData(histogramHarmonic.getData());
		        dataLogger.closeFile();
		        /*
		         * 
		         */
		
		        double wHarmonic = dataAverageHarmonic.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
		        double wTarget   = dataAverageTarget.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);        
		      
		        /*
		         * Qharmonic = < e0 / [sqrt(e1^2 + alpha^2 * e0^2)]>umbrella
		         *  Qtarget  = < e1 / [sqrt(e1^2 + alpha^2 * e0^2)]>umbrella
		         * 
		         */
		        double Qharmonic = dataAverageSamplingHarmonic.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
		        double Qtarget = dataAverageSamplingTarget.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
		
		        /*
		         * deltaFE_harmonic: beta*(FE_harmonic - FE_umbrella) = - ln(Qharmonic)
		         *  deltaFE_target : beta*(FE_target - FE_umbrella) = - ln(Qtarget)
		         */
		        double deltaFE_harmonic = - Math.log(Qharmonic);
		        double deltaFE_target = - Math.log(Qtarget);
		        double deltaFE = temp*deltaFE_target - temp*deltaFE_harmonic;
		        
		        double sHarmonic = wHarmonic - deltaFE_harmonic;
		        double sTarget = wTarget - deltaFE_target;
	
		        double reweightedwHarmonic = meterWorkUmbrellaHarmonic.getDataReweighted();
		        double reweightedwTarget = meterWorkUmbrellaTarget.getDataReweighted();

		        double pi_harmonic = Math.sqrt((sHarmonic/sTarget)*Math.log((0.5/Math.PI)*(steps-1)*(steps-1)))-Math.sqrt(2*sHarmonic);
		        double pi_target = Math.sqrt((sTarget/sHarmonic)*Math.log((0.5/Math.PI)*(steps-1)*(steps-1)))-Math.sqrt(2*sTarget);
	
		        
		        try {
		        	fileWriterUmb.write( idStep + " " + deltaFE_harmonic + " " + deltaFE_target + " " + deltaFE + " "
		        			                          + wHarmonic + " " + wTarget + " "
		        			                          + reweightedwHarmonic + " " + reweightedwTarget + " "
		        			                          + sHarmonic + " " + sTarget + " "
		        			                          + pi_harmonic + " " + pi_target + "\n");
		        } catch (IOException e){
		        	
		        }
		        
        	}
        };
        
        sim.integrator.addIntervalAction(outputAction);
        sim.integrator.setActionInterval(outputAction, 10000);
        
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.getController().actionPerformed();
        
        try{
	        fileWriterUmb.close();
        } catch (IOException e){
        	
        }

		double[][] omega2 = sim.normalModes.getOmegaSquared(sim.box);
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
		AHarmonic -= Math.log(Math.pow(2.0, basisSize*D*(totalCells - fac)/2.0) / Math.pow(totalCells,0.5*D));
		System.out.println("Harmonic-reference free energy: "+AHarmonic*temperature);
        
        double wHarmonic = dataAverageHarmonic.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
        double wTarget   = dataAverageTarget.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);        
        
        double eHarmonic = dataAverageHarmonic.getData().getValue(AccumulatorAverage.StatType.ERROR.index);
        double eTarget   = dataAverageTarget.getData().getValue(AccumulatorAverage.StatType.ERROR.index);
        /*
         * Qharmonic = < e0 / [sqrt(e1^2 + alpha^2 * e0^2)]>umbrella
         *  Qtarget  = < e1 / [sqrt(e1^2 + alpha^2 * e0^2)]>umbrella
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
        double deltaFE = temp*deltaFE_target - temp*deltaFE_harmonic;
        
        double sHarmonic = wHarmonic - deltaFE_harmonic;
        double sTarget = wTarget - deltaFE_target;
        
        double er_sHarmonic = Math.sqrt(eHarmonic*eHarmonic + (eQharmonic/Qharmonic)*(eQharmonic/Qharmonic));
        double er_sTarget   = Math.sqrt(eTarget*eTarget + (eQtarget/Qtarget)*(eQtarget/Qtarget));
        
        double reweightedwHarmonic = meterWorkUmbrellaHarmonic.getDataReweighted();
        double reweightedwTarget = meterWorkUmbrellaTarget.getDataReweighted();
        
        System.out.println("\nfree energy difference (umbrella --> harmonic): "+ (temp*deltaFE_harmonic) + "; error: " + eQharmonic/Qharmonic);
        System.out.println("free energy difference (umbrella --> target): "+ (temp*deltaFE_target) + "; error: " + eQtarget/Qtarget);
        System.out.println("free energy difference (harmonic --> target): " + deltaFE 
	           	+ ", error: " + Math.sqrt( (eQharmonic/Qharmonic)*(eQharmonic/Qharmonic) + (eQtarget/Qtarget)*(eQtarget/Qtarget) ));
        System.out.println("target free energy: " + (temperature*AHarmonic+ deltaFE)+ 
        		" ,error: " + Math.sqrt( (eQharmonic/Qharmonic)*(eQharmonic/Qharmonic) + (eQtarget/Qtarget)*(eQtarget/Qtarget) ));
        System.out.println("target free energy per particle: " + (temperature*AHarmonic+ deltaFE)/numAtoms);
        
        
        System.out.println("\nwUmbrellaHarmonic: "+ wHarmonic + ", error: "+ eHarmonic);
        System.out.println(" wUmbrellaTarget : "  + wTarget   + ", error: "+ eTarget);
        
        System.out.println("\nwHarmonicUmbrella: "+ reweightedwHarmonic );
        System.out.println(" wTargetUmbrella : "  + reweightedwTarget   );
        
        System.out.println(" beta*deltaFE (harmonic): " + deltaFE_harmonic);
        System.out.println(" beta*deltaFE (target): " + deltaFE_target);
        
        System.out.println("\nUmbrella (perturbed into Harmonic) entropy: " + sHarmonic + ", error: "+ er_sHarmonic);
        System.out.println("Umbrella (perturbed into Target) entropy: " + sTarget + ", error: "+ er_sTarget);
        
        double pi_harmonic = Math.sqrt((sHarmonic/sTarget)*Math.log((0.5/Math.PI)*(steps-1)*(steps-1)))-Math.sqrt(2*sHarmonic);
        double pi_target = Math.sqrt((sTarget/sHarmonic)*Math.log((0.5/Math.PI)*(steps-1)*(steps-1)))-Math.sqrt(2*sTarget);
        
        System.out.println("PI Harmonic: " + pi_harmonic);
        System.out.println("PI Target: " + pi_target);
        
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
    public double refPref;
    
    public static class SimBennetParam extends ParameterBase {
    	public int numMolecules = 32;
    	public double density = 1256;
    	public int exponentN = 12;
    	public int D = 3;
    	public long numSteps = 2000000;
    	public double harmonicFudge =1;
    	public String filename = "CB_FCC_n12_T01";
    	public double temperature = 0.1;
    }

}