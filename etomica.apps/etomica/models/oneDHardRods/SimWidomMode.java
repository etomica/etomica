package etomica.models.oneDHardRods;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomList;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.atom.Atom;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPump;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.listener.IntegratorListenerAction;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinitionLeaf;
import etomica.normalmode.MCMoveAtomCoupled;
import etomica.normalmode.NormalModes;
import etomica.normalmode.NormalModes1DHR;
import etomica.normalmode.P2XOrder;
import etomica.normalmode.WaveVectorFactory;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2;
import etomica.potential.Potential2HardSpherical;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.RandomNumberGenerator;
import etomica.util.ReadParameters;

/**
 * MD simulation of hard spheres in 1D or 3D with tabulation of the
 * collective-coordinate S-matrix. No graphic display of simulation.
 */
public class SimWidomMode extends Simulation {

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "SimDegreeFreedom1DHR";
    public Primitive primitive;
    int[] nCells;
    NormalModes nm;
    public IntegratorMC integrator;
    public BasisMonatomic basis;
    public ActivityIntegrate activityIntegrate;
    
    public IBox box;
    public Boundary bdry;
    public CoordinateDefinition coordinateDefinition;
    WaveVectorFactory waveVectorFactory;
    MCMoveAtomCoupled mcMoveAtom;
    MCMoveChangeSingleMode mcMoveMode;
    int harmonicWV;
    MeterWidomMode[] widomMeter;
    AccumulatorAverage[] accumulators;

    public SimWidomMode(Space _space, int numAtoms, double density, int blocksize) {
        super(_space, true);
        
//        long seed = 3;
//        System.out.println("Seed explicitly set to " + seed);
//        IRandom rand = new RandomNumberGenerator(seed);
//        this.setRandom(rand);
        
        PotentialMasterList potentialMaster = new PotentialMasterList(this, space);

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);
        basis = new BasisMonatomic(space);
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
       
        Potential2 potential = new P2HardSphere(space, 1.0, true);
        potential = new P2XOrder(space, (Potential2HardSpherical)potential);
        potential.setBox(box);
        potentialMaster.addPotential(potential, new IAtomType[] {species.getLeafType(), species.getLeafType()});

        primitive = new PrimitiveCubic(space, 1.0/density);
        bdry = new BoundaryRectangularPeriodic(space, numAtoms/density);
        nCells = new int[]{numAtoms};
        box.setBoundary(bdry);
        
        coordinateDefinition = new CoordinateDefinitionLeaf(this, box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(nCells);
        
        double neighborRange = 1.01/density;
        potentialMaster.setRange(neighborRange);
        //find neighbors now.  Don't hook up NeighborListManager since the
        //  neighbors won't change
        potentialMaster.getNeighborManager(box).reset();
        
        integrator = new IntegratorMC(this, potentialMaster);
        integrator.setBox(box);
        
        nm = new NormalModes1DHR(space.D());
        nm.setHarmonicFudge(1.0);
        nm.setTemperature(1.0);
        nm.getOmegaSquared(box);
        
        waveVectorFactory = nm.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);
        
        mcMoveAtom = new MCMoveAtomCoupled(potentialMaster, random, space);
        mcMoveAtom.setPotential(potential);
        mcMoveAtom.setBox(box);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        mcMoveAtom.setStepSizeMin(0.001);
        mcMoveAtom.setStepSize(0.01);
        
        mcMoveMode = new MCMoveChangeSingleMode(potentialMaster, random);
        mcMoveMode.setBox(box);
        integrator.getMoveManager().addMCMove(mcMoveMode);
        mcMoveMode.setCoordinateDefinition(coordinateDefinition);
        mcMoveMode.setEigenVectors(nm.getEigenvectors(box));
        mcMoveMode.setOmegaSquared(nm.getOmegaSquared(box));
        mcMoveMode.setWaveVectorCoefficients(nm.getWaveVectorFactory().getCoefficients());
        mcMoveMode.setWaveVectors(nm.getWaveVectorFactory().getWaveVectors());
        
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        int coordNum = nm.getWaveVectorFactory().getWaveVectors().length*coordinateDim;
        
        widomMeter = new MeterWidomMode[coordNum];
        accumulators = new AccumulatorAverageFixed[coordNum];
        DataPump pump;
        IntegratorListenerAction pumpListener;
        for(int i = 0; i < coordNum; i++){
            String name = new String("widom Meter for mode " + i);
            widomMeter[i] = new MeterWidomMode(name, potentialMaster, 
                    coordinateDefinition, box, i);
            widomMeter[i].setEigenVectors(nm.getEigenvectors(box));
            widomMeter[i].setOmegaSquared(nm.getOmegaSquared(box));
            widomMeter[i].setWaveVectorCoefficients(nm.getWaveVectorFactory().getCoefficients());
            widomMeter[i].setWaveVectors(nm.getWaveVectorFactory().getWaveVectors());
            
            accumulators[i] = new AccumulatorAverageFixed(blocksize);
            
            pump = new DataPump(widomMeter[i], accumulators[i]);
            pumpListener = new IntegratorListenerAction(pump);
            pumpListener.setInterval(blocksize);
            integrator.getEventManager().addListener(pumpListener);
        }
        
        activityIntegrate = new ActivityIntegrate(integrator, 0, true);
        getController().addAction(activityIntegrate);
        
        
//        IAtomList leaflist = box.getLeafList();
//        double[] locations = new double[numAtoms];
//        System.out.println("at sim start:");
//        for(int i = 0; i < numAtoms; i++){
//            //one d is assumed here.
//            locations[i] = ( ((Atom)leaflist.getAtom(i)).getPosition().x(0) );
//        }
//        for(int i = 0; i < numAtoms; i++){
//            System.out.println(i + "  " + locations[i]);
//        }
    }

    private void setHarmonicWV(int hwv){
        harmonicWV = hwv;
        mcMoveMode.setHarmonicWV(hwv);
    }
    
    /**
     * @param args
     */
    public static void main(String[] args) {

        SimParam params = new SimParam();
        String inputFilename = null;
        if(args.length > 0) {
            inputFilename = args[0];
        }
        if(inputFilename != null){
            ReadParameters readParameters = new ReadParameters(inputFilename, params);
            readParameters.readParameters();
            inputFilename = params.inputfilename;
        }
        
        int nA = params.numAtoms;
        double density = params.density;
        int D = params.D;
        double harmonicFudge = params.harmonicFudge;
        String filename = params.filename;
        if(filename.length() == 0){
            filename = "1DHR";
        }
        double temperature = params.temperature;
        int comparedWV = params.comparedWV;
        int nSteps = params.numSteps;
        int bs = params.blockSize;
        
        System.out.println("Running "
                + (D == 1 ? "1D" : (D == 3 ? "FCC" : "2D hexagonal"))
                + " hard sphere simulation");
        System.out.println(nA + " atoms at density " + density);
        System.out.println(nSteps + " steps, " + bs + " blocksize");
        System.out.println("input data from " + inputFilename);
        System.out.println("output data to " + filename);

        SimWidomMode sim = new SimWidomMode(Space.getInstance(D), nA, density, bs);
        
        // start simulation
        sim.activityIntegrate.setMaxSteps(nSteps/10);
        sim.setHarmonicWV(comparedWV);
        sim.getController().actionPerformed();
        System.out.println("equilibration finished");
        sim.getController().reset();
       
        sim.activityIntegrate.setMaxSteps(nSteps);
        sim.getController().actionPerformed();
        
        //After processing...
        int cd = sim.nm.getWaveVectorFactory().getWaveVectors().length * 
                sim.coordinateDefinition.getCoordinateDim();
        double[] results = new double[cd];
        DataGroup group;
        for(int i = 0; i < cd; i++){
            group = (DataGroup)sim.accumulators[i].getData();
            results[i] = ((DataDouble)group.getData(StatType.AVERAGE.index)).x;
        }
        for(int i = 0; i < cd; i++){
            System.out.println(i + "  " + results[i]);
        }
        
//        IAtomList leaflist = sim.box.getLeafList();
//        double[] locations = new double[nA];
//        System.out.println("final:");
//        for(int i = 0; i < nA; i++){
//            //one d is assumed here.
//            locations[i] = ( ((Atom)leaflist.getAtom(i)).getPosition().x(0) );
//            System.out.println(i + "  " + locations[i]);
//        }
        
        System.out.println("Fini.");
    }
    
    public static class SimParam extends ParameterBase {
        public int numAtoms = 32;
        public double density = 0.50;
        public int D = 1;
        public double harmonicFudge = 1.0;
        public String filename = "HR1D_";
        public String inputfilename = "input";
        public double temperature = 1.0;
        public int comparedWV = numAtoms/2;
        
        public int blockSize = 100;
        public int numSteps = 10000000;
    }

}