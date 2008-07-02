package etomica.normalmode;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistogram;
import etomica.data.DataLogger;
import etomica.data.DataPump;
import etomica.data.DataSource;
import etomica.data.DataTableWriter;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncated;
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
 * Simulation to sample Umbrella Overlap Sampling
 */
public class SimUmbrella extends Simulation {

	private static final String APP_NAME = "Sim Bennet's";

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
        normalModes.setTemperature(temperature);
        coordinateDefinition.initializeCoordinates(nCells);
        
        Potential2SoftSpherical potential = new P2SoftSphere(space, 1.0, 1.0, exponent);
        double truncationRadius = boundary.getDimensions().x(0) * 0.495;
        P2SoftSphericalTruncated pTruncated = new P2SoftSphericalTruncated(space, potential, truncationRadius);
        IAtomTypeLeaf sphereType = species.getLeafType();
        potentialMasterMonatomic.addPotential(pTruncated, new IAtomTypeLeaf[] { sphereType, sphereType });
        
        integrator.setBox(box);
        
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
        SimUmbrella sim = new SimUmbrella(Space.getInstance(D), numAtoms, density, temperature, filename, exponentN);
        
        DataSource[] workMeters = new DataSource[2];
        DataSource[] samplingMeters = new DataSource[2];
        
      // Umbrella Sampling ---> Harmonic
        MeterWorkUmbrellaHarmonic meterWorkUmbrellaHarmonic = new MeterWorkUmbrellaHarmonic(sim.integrator, sim.meterHarmonicEnergy);
        meterWorkUmbrellaHarmonic.setRefPref(sim.refPref);
        meterWorkUmbrellaHarmonic.setLatticeEnergy(sim.latticeEnergy);
        workMeters[0] = meterWorkUmbrellaHarmonic;
        
        AccumulatorAverageFixed dataAverageHarmonic = new AccumulatorAverageFixed();
        DataPump pumpHarmonic = new DataPump(workMeters[0], dataAverageHarmonic);
        sim.integrator.addIntervalAction(pumpHarmonic);
        sim.integrator.setActionInterval(pumpHarmonic, 1);
        
        //Histogram Harmonic
        AccumulatorHistogram histogramHarmonic = new AccumulatorHistogram(new HistogramExpanding(1, new DoubleRange(0,1)));
        DataPump pumpHistogramHarmonic = new DataPump(workMeters[0], histogramHarmonic);
        sim.integrator.addIntervalAction(pumpHistogramHarmonic);
        sim.integrator.setActionInterval(pumpHistogramHarmonic, 1);
        
      // Umbrella Sampling ---> Target
        MeterWorkUmbrellaTarget meterWorkUmbrellaTarget = new MeterWorkUmbrellaTarget(sim.integrator, sim.meterHarmonicEnergy);
        meterWorkUmbrellaTarget.setRefPref(sim.refPref);
        meterWorkUmbrellaTarget.setLatticeEnergy(sim.latticeEnergy);
        workMeters[1] = meterWorkUmbrellaTarget;
        
        AccumulatorAverageFixed dataAverageTarget = new AccumulatorAverageFixed();
        DataPump pumpTarget = new DataPump(workMeters[1], dataAverageTarget);
        sim.integrator.addIntervalAction(pumpTarget);
        sim.integrator.setActionInterval(pumpTarget, numAtoms*2);
        
        //Histogram Target
        AccumulatorHistogram histogramTarget = new AccumulatorHistogram(new HistogramExpanding(1,new DoubleRange(0,1)));
        DataPump pumpHistogramTarget = new DataPump(workMeters[1], histogramTarget);
        sim.integrator.addIntervalAction(pumpHistogramTarget);
        sim.integrator.setActionInterval(pumpHistogramTarget, numAtoms*2); //the interval is based on the system size
        

        /*
         * 
         */
        
        // Harmonic Sampling
        
        MeterSamplingHarmonic meterSamplingHarmonic = new MeterSamplingHarmonic(sim.integrator, sim.meterHarmonicEnergy);
        meterSamplingHarmonic.setRefPref(sim.refPref);
        meterSamplingHarmonic.setLatticeEnergy(sim.latticeEnergy);
        samplingMeters[0] = meterSamplingHarmonic;
        
        AccumulatorAverageFixed dataAverageSamplingHarmonic = new AccumulatorAverageFixed();
        DataPump pumpSamplingHarmonic = new DataPump(samplingMeters[0], dataAverageSamplingHarmonic);
        sim.integrator.addIntervalAction(pumpSamplingHarmonic);
        sim.integrator.setActionInterval(pumpSamplingHarmonic, 1);
        
        // Target Sampling
        
        MeterSamplingTarget meterSamplingTarget = new MeterSamplingTarget(sim.integrator, sim.meterHarmonicEnergy);
        meterSamplingTarget.setRefPref(sim.refPref);
        meterSamplingTarget.setLatticeEnergy(sim.latticeEnergy);
        samplingMeters[1] = meterSamplingTarget;
        
        AccumulatorAverageFixed dataAverageSamplingTarget = new AccumulatorAverageFixed();
        DataPump pumpSamplingTarget = new DataPump(samplingMeters[1], dataAverageSamplingTarget);
        sim.integrator.addIntervalAction(pumpSamplingTarget);
        sim.integrator.setActionInterval(pumpSamplingTarget, 1);
        
        
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.getController().actionPerformed();
        
        /*
         * Histogram
         */
        //Target
		DataLogger dataLogger = new DataLogger();
		DataTableWriter dataTableWriter = new DataTableWriter();
		dataLogger.setFileName(filename + "_hist_UmbreTarg");
		dataLogger.setDataSink(dataTableWriter);
		dataTableWriter.setIncludeHeader(false);
		dataLogger.putDataInfo(histogramTarget.getDataInfo());
		
		dataLogger.setWriteInterval(1);
		dataLogger.putData(histogramTarget.getData());
		dataLogger.closeFile();
        
		//Harmonic
        dataLogger.setFileName(filename + "_hist_UmbreHarm");
        dataLogger.putData(histogramHarmonic.getData());
        dataLogger.closeFile();
        /*
         * 
         */
        
        
        
        
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
         * deltaFE_harmonic: beta*(U_harmonic - U_umbrella) = - ln(Qharmonic)
         *  deltaFE_target : beta*(U_target - U_umbrella) = - ln(Qtarget)
         */
        double deltaFE_harmonic = - Math.log(Qharmonic);
        double deltaFE_target = - Math.log(Qtarget);
        
        double sHarmonic = wHarmonic - deltaFE_harmonic;
        double sTarget = wTarget - deltaFE_target;
        
        double er_sHarmonic = Math.sqrt(eHarmonic*eHarmonic + (eQharmonic/Qharmonic)*(eQharmonic/Qharmonic));
        double er_sTarget   = Math.sqrt(eTarget*eTarget + (eQtarget/Qtarget)*(eQtarget/Qtarget));
        
        System.out.println("\nfree energy difference (harmonic): "+ (temperature*deltaFE_harmonic) + "; error: " + eQharmonic/Qharmonic);
        System.out.println("free energy difference (target): "+ (temperature*deltaFE_target) + "; error: " + eQtarget/Qtarget);
        
        System.out.println("\nwUmbrellaHarmonic: "+ wHarmonic + ", error: "+ eHarmonic);
        System.out.println(" wUmbrellaTarget : "  + wTarget   + ", error: "+ eTarget);
        System.out.println(" deltaFE (harmonic): " + deltaFE_harmonic);
        System.out.println(" deltaFE (target): " + deltaFE_target);
        
        System.out.println("\nUmbrella (perturbed into Harmonic) entropy: " + sHarmonic + ", error: "+ er_sHarmonic);
        System.out.println("Umbrella (perturbed into Target) entropy: " + sTarget + ", error: "+ er_sTarget);
        
        //double pi_harmonic = Math.sqrt((sHarmonic/sTarget)*Math.log((0.5/Math.PI)*numSteps*numSteps))-Math.sqrt(2*sHarmonic);
        //double pi_target = Math.sqrt((sTarget/sHarmonic)*Math.log((0.5/Math.PI)*numSteps*numSteps))-Math.sqrt(2*sTarget);
        
        //System.out.println("PI Harmonic: " + pi_harmonic);
        //System.out.println("PI Target: " + pi_target);
        
        
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
    	public long numSteps = 100000;
    	public double harmonicFudge =1;
    	public String filename = "FCC_SoftSphere_n12_T02";
    	public double temperature = 0.2;
    }

}