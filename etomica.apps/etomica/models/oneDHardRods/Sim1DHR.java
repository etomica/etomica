package etomica.models.oneDHardRods;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.data.DataPump;
import etomica.data.IEtomicaDataSource;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinitionLeaf;
import etomica.normalmode.NormalModes1DHR;
import etomica.normalmode.P2XOrder;
import etomica.normalmode.SimTarget;
import etomica.normalmode.WaveVectorFactory;
import etomica.potential.P1HardPeriodic;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2;
import etomica.potential.Potential2HardSpherical;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;

/**
 * Simulation that changes one normal mode to a Gaussian defined position
 * 
 * @author cribbin
 */
public class Sim1DHR extends Simulation {

    private static final long serialVersionUID = 1L;
    public DataSourceVirialOverlap dsvo;
    public Boundary boundaryTarget, boundaryOriginal;
    public double refPref;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    public DataPump[] accumulatorPumps;
    public IEtomicaDataSource[] meters;
	   
	Integrator integratorOverlap;
	IntegratorBox[] integrators;
	ActivityIntegrate activityIntegrate;
	IBox boxTarget, boxOriginal;
	Boundary boundary;
	CoordinateDefinition coordinateDefinition;
	Primitive primitive;
	Basis basis;
	int[] nCells;
	SpeciesSpheresMono species;
	NormalModes1DHR nm;
	
	private static final String APP_NAME = "Sim 1DHR";

	public Sim1DHR(Space _space, int numAtoms, double density, double 
			temperature, String filename, double harmonicFudge){
        super(_space, true);
     
        integrators = new IntegratorBox[2];
        accumulatorPumps = new DataPump[2];
        meters = new IEtomicaDataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];
        
        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);
        
        //TARGET - 1DHR system
        PotentialMasterList potentialMasterTarget = new PotentialMasterList(this, space);
        boxTarget = new Box(space);
        boxTarget.setNMolecules(species, numAtoms);
        
        IntegratorHard integratorTarget = new IntegratorHard(this, potentialMasterTarget, space);
        integratorTarget.setTimeStep(4);
        integratorTarget.setTemperature(temperature);
        integratorTarget.setIsothermal(false);
        integrators[1] = integratorTarget;
        
        Potential2 p2 = new P2HardSphere(space, 1.0, true);
        p2 = new P2XOrder(space, (Potential2HardSpherical)p2);
        potentialMasterTarget.addPotential(p2, new IAtomTypeLeaf[]
                {species.getLeafType(), species.getLeafType()});
        
        primitive = new PrimitiveCubic(space, 1.0/density);
        boundaryTarget = new BoundaryRectangularPeriodic(space, numAtoms/density);
        integratorTarget.setNullPotential(new P1HardPeriodic(space), species.getLeafType());
        nCells = new int[]{numAtoms};
        boxTarget.setBoundary(boundaryTarget);
        
        CoordinateDefinitionLeaf coordinateDefinitionTarget = new 
                CoordinateDefinitionLeaf(this, boxTarget, primitive, space);
        coordinateDefinitionTarget.initializeCoordinates(nCells);
        
        double neighborRange = 1.01/density;
        potentialMasterTarget.setRange(neighborRange);
        //find neighbors now.  Don't hook up NeighborListManager since the
        //  neighbors won't change
        potentialMasterTarget.getNeighborManager(boxTarget).reset();
        integratorTarget.setBox(boxTarget);
        
        
        //ORIGINAL
        boundaryOriginal = new BoundaryRectangularPeriodic(space, numAtoms/density);
        boxOriginal = new Box(boundaryOriginal, space);
        addBox(boxOriginal);
        boxOriginal.setNMolecules(species,numAtoms);
        IntegratorMC integratorOriginal = new IntegratorMC(null, random, temperature);
        integratorOriginal.setBox(boxOriginal);
        
        MCMoveChangeMode convert = new MCMoveChangeMode(potentialMasterTarget, random);
        integratorOriginal.getMoveManager().addMCMove(convert);
        integrators[0] = integratorOriginal;
        
        CoordinateDefinitionLeaf coordinateDefinitionOriginal = new 
                CoordinateDefinitionLeaf(this, boxOriginal, primitive, space);
        
        nm = new NormalModes1DHR(space.D());
        nm.setHarmonicFudge(harmonicFudge);
        nm.setTemperature(temperature);
        
        WaveVectorFactory waveVectorFactory = nm.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(boxOriginal);
        
        //nan WAVE VECTOR FACTORY STUFF GOES HERE
        
        convert.setCoordinateDefinition(coordinateDefinitionOriginal);
//        convert.setTemperature(temperature);
        convert.setBox((IBox)boxOriginal);
        
        integratorOriginal.setBox(boxOriginal);
        potentialMasterTarget.getNeighborManager(boxOriginal).reset();
        
        
        //OVERLAP   
        
        
        
        
        
        
        
       

        
        
        nm = new NormalModes1DHR(1);
        
        
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
 
		/*
		 * This whole setup defines a set of default parameters
		 * in the inner class Sim1DHRParams.  These parameters can be changed
		 * individually in an appropriately named file, without affecting
		 * the values of the other parameters.  The order of definition in the
		 * file is irrelevant.
		 * 
		 */
		Sim1DHRParams params = new Sim1DHRParams();
		String inputFilename = null;
		if(args.length > 0){
			inputFilename = args[0];
		}
        if(inputFilename != null){
        	ReadParameters readParameters = new ReadParameters(inputFilename, params);
        	readParameters.readParameters();
        }
        double density = params.density;
        double numSteps = params.numSteps;
        int numAtoms = params.numAtoms;
        double harmonicFudge = params.harmonicFudge;
        double temperature = params.temperature;
        int D = params.D;
        String filename = params.filename;
        if(filename.length() ==0){
        	filename = "normal_modes_1DHR _" + numAtoms;
        }
        String refFileName = args.length>0 ? filename+"_ref" : null;
        
        System.out.println("Running 1D hard rod simulation");
        System.out.println(numAtoms+" atoms at density "+density);
        System.out.println("harmonic fudge: "+harmonicFudge);
        System.out.println((numSteps/1000)+ " total steps of 1000");
        System.out.println("output data to "+filename);
        
        
        
        
        //instantiate simulation
        SimTarget sim = new SimTarget(Space.getInstance(D), numAtoms, density);
        
        
        
        
        
        
        
        System.out.println("Fini.");
	}

	
    /**
     * Inner class for parameters understood by the Sim1DHR constructor
     */
    public static class Sim1DHRParams extends ParameterBase {
        public int numAtoms = 32;
        public double density = 0.5;
        public int D = 1;
        public long numSteps = 100000;
        public double harmonicFudge = 1.0;
        public String filename = "HR1D_";
        public double temperature = 1.0;
    }
}
