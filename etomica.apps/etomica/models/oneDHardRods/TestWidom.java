package etomica.models.oneDHardRods;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistogram;
import etomica.data.DataPump;
import etomica.data.DataSplitter;
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
import etomica.util.DoubleRange;
import etomica.util.Histogram;
import etomica.util.HistogramSimple;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;

/**
 * MD simulation of hard spheres in 1D or 3D with tabulation of the
 * collective-coordinate S-matrix. No graphic display of simulation.
 */
public class TestWidom extends Simulation {

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
    MCMoveChangeSingleWVLoop mcMoveMode;
    int harmonicWV, pickedWV;
    MeterWidomModeReal realMeter;
    MeterWidomModeImaginary imagMeter;
    MeterNormalModeCoordinate mnm;
    AccumulatorAverage realAccumulator, imagAccumulator;
    AccumulatorHistogram[] hists;

    public TestWidom(Space _space, int numAtoms, double density, int blocksize, int nbs, int pwv) {
        super(_space, true);
        
        pickedWV = pwv;
        
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
        
        nm = new NormalModes1DHR(bdry, numAtoms);
        nm.setHarmonicFudge(1.0);
        nm.setTemperature(1.0);
        nm.getOmegaSquared();
        
        waveVectorFactory = nm.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);
        
        mcMoveAtom = new MCMoveAtomCoupled(potentialMaster, random, space);
        mcMoveAtom.setPotential(potential);
        mcMoveAtom.setBox(box);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        mcMoveAtom.setStepSizeMin(0.001);
        mcMoveAtom.setStepSize(0.01);
        
        mcMoveMode = new MCMoveChangeSingleWVLoop(potentialMaster, random);
        mcMoveMode.setBox(box);
        integrator.getMoveManager().addMCMove(mcMoveMode);
        mcMoveMode.setCoordinateDefinition(coordinateDefinition);
        mcMoveMode.setEigenVectors(nm.getEigenvectors());
        mcMoveMode.setOmegaSquared(nm.getOmegaSquared());
        mcMoveMode.setWaveVectorCoefficients(nm.getWaveVectorFactory().getCoefficients());
        mcMoveMode.setWaveVectors(nm.getWaveVectorFactory().getWaveVectors());
        
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        int coordNum = nm.getWaveVectorFactory().getWaveVectors().length*coordinateDim * 2;
        
        String name = new String("widom Meter for real mode " + pickedWV);
        realMeter = new MeterWidomModeReal(name, potentialMaster, 
                coordinateDefinition, box, pickedWV);
        realMeter.setEigenVectors(nm.getEigenvectors());
        realMeter.setOmegaSquared(nm.getOmegaSquared());
        realMeter.setWaveVectorCoefficients(nm.getWaveVectorFactory().getCoefficients());
        realMeter.setWaveVectors(nm.getWaveVectorFactory().getWaveVectors());
        
        realAccumulator = new AccumulatorAverageFixed(blocksize);
        DataPump pump = new DataPump(realMeter, realAccumulator);
        IntegratorListenerAction pumpListener;
        pumpListener = new IntegratorListenerAction(pump);
        pumpListener.setInterval(blocksize);
        integrator.getEventManager().addListener(pumpListener);
        
        
        name = new String("widom Meter for imag mode " + pickedWV);
        imagMeter = new MeterWidomModeImaginary(name, potentialMaster, 
                coordinateDefinition, box, pickedWV);
        imagMeter.setEigenVectors(nm.getEigenvectors());
        imagMeter.setOmegaSquared(nm.getOmegaSquared());
        imagMeter.setWaveVectorCoefficients(nm.getWaveVectorFactory().getCoefficients());
        imagMeter.setWaveVectors(nm.getWaveVectorFactory().getWaveVectors());
        
        imagAccumulator = new AccumulatorAverageFixed(blocksize);
        pump = new DataPump(imagMeter, imagAccumulator);
        pumpListener = new IntegratorListenerAction(pump);
        pumpListener.setInterval(blocksize);
        integrator.getEventManager().addListener(pumpListener);
        
        mnm = new MeterNormalModeCoordinate(coordinateDefinition, nm.getWaveVectorFactory().getWaveVectors());
        mnm.setEigenVectors(nm.getEigenvectors());
        mnm.setOmegaSquared(nm.getOmegaSquared());
        
        hists = new AccumulatorHistogram[coordNum];
        DataSplitter splitter = new DataSplitter();
        pump = new DataPump(mnm, splitter);
        
        DoubleRange range = new DoubleRange(-1.0, 1.0);
        Histogram template;
        for(int i = 0; i < coordNum; i++){
            template = new HistogramSimple(nbs, range);
            hists[i] = new AccumulatorHistogram(template, nbs);
            splitter.setDataSink(i, hists[i]);
        }
        
        pumpListener = new IntegratorListenerAction(pump);
        pumpListener.setInterval(blocksize);
        integrator.getEventManager().addListener(pumpListener);
        
        activityIntegrate = new ActivityIntegrate(integrator, 0, true);
        getController().addAction(activityIntegrate);
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
        int nbins = params.nBins;
        int picked = params.pickedWV;
        
        System.out.println("Running "
                + (D == 1 ? "1D" : (D == 3 ? "FCC" : "2D hexagonal"))
                + " hard sphere simulation");
        System.out.println(nA + " atoms at density " + density);
        System.out.println(nSteps + " steps, " + bs + " blocksize");
        System.out.println("input data from " + inputFilename);
        System.out.println("output data to " + filename);
        System.out.println("picked WV " + picked);

        TestWidom sim = new TestWidom(Space.getInstance(D), nA, density, bs, nbins, picked);
        
        // start simulation
        sim.activityIntegrate.setMaxSteps(nSteps/10);
        sim.setHarmonicWV(comparedWV);
        sim.getController().actionPerformed();
        System.out.println("equilibration finished");
        

        int accumulatorLength = sim.hists.length;
        for(int i = 0; i < accumulatorLength; i++){
            sim.hists[i].reset();
        }
        sim.getController().reset();
       
        sim.activityIntegrate.setMaxSteps(nSteps);
        sim.getController().actionPerformed();
        
        //After processing...
        DataGroup realgroup = (DataGroup)sim.realAccumulator.getData();
        double realresult = ((DataDouble)realgroup.getData(StatType.AVERAGE.index)).x;
        System.out.println("real results " + realresult);
         
        DataGroup imaggroup = (DataGroup)sim.imagAccumulator.getData();
        double imagresult = ((DataDouble)imaggroup.getData(StatType.AVERAGE.index)).x;
        System.out.println("imag results " + imagresult);
         
        /* 
         * This loop creates a new write class for each histogram from each 
         * AccumulatorHistogram, changes the filename for the histogram output,
         * connects the write class with this histogram, and 
         * writes out the results to the file.
         */
        String outputfileName = new String("hist_" + sim.pickedWV);
        WriteHistograms wh = new WriteHistograms(outputfileName);
        wh.setHistogram(sim.hists[sim.pickedWV].getHistograms());
        wh.actionPerformed();
        System.out.println("hist count " + sim.hists[sim.pickedWV].getHistograms().getCount());
        
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
        
        public int blockSize = 1000;
        public int numSteps = 10000000;
        public int nBins = 200;
        public int pickedWV = 12;
    }

}