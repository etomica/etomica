package etomica.models.oneDHardRods;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.data.AccumulatorHistogram;
import etomica.data.DataPump;
import etomica.data.DataSplitter;
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
import etomica.util.HistogramExpanding;

/**
 * MD simulation of hard spheres in 1D or 3D with tabulation of the
 * collective-coordinate S-matrix. No graphic display of simulation.
 */
public class SimDegreeFreedom extends Simulation {

    public SimDegreeFreedom(Space _space, int numAtoms, double density) {
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
//        AtomTypeSphere sphereType = (AtomTypeSphere)species.getLeafType();
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
        
        mcmove = new MCMoveAtomCoupled(potentialMaster, random, space);
        mcmove.setPotential(potential);
        mcmove.setBox(box);
        integrator.getMoveManager().addMCMove(mcmove);
        mcmove.setStepSizeMin(0.001);
        mcmove.setStepSize(0.01);
        
        meternmc = new MeterNormalModeCoordinate(coordinateDefinition, nm.getWaveVectorFactory().getWaveVectors());
        meternmc.setEigenVectors(nm.getEigenvectors(box));
        meternmc.setOmegaSquared(nm.getOmegaSquared(box));
        
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        int coordNum = nm.getWaveVectorFactory().getWaveVectors().length*coordinateDim*2;
        hists = new AccumulatorHistogram[coordNum];
        DataSplitter splitter = new DataSplitter();
        DataPump pumpFromMeter = new DataPump(meternmc, splitter);
        
        DoubleRange range = new DoubleRange(-1.0, 1.0);
        int nBins = 200;
        HistogramExpanding template; 
        for(int i = 0; i < coordNum; i++){
            template = new HistogramExpanding(nBins, range);
            hists[i] = new AccumulatorHistogram(template, nBins);
            splitter.setDataSink(i, hists[i]);
        }
        
        IntegratorListenerAction pumpFromMeterListener = new IntegratorListenerAction(pumpFromMeter);
        pumpFromMeterListener.setInterval(1000);
        integrator.getEventManager().addListener(pumpFromMeterListener);
        
        activityIntegrate = new ActivityIntegrate(integrator, 0, true);
        getController().addAction(activityIntegrate);
        
        
    }

    /**
     * @param args
     */
    public static void main(String[] args) {

        // defaults
        int D = 1;
        int nA = 32;
        double density = 1.3;
        if (D == 1) {
            nA = 32;
            density = 0.5;
        }

        int nSteps = 10000000;

        // parse arguments
        if (args.length > 1) {
            density = Double.parseDouble(args[1]);
        }
        if (args.length > 2) {
//            simTime = Double.parseDouble(args[2]);
        }
        if (args.length > 3) {
            nA = Integer.parseInt(args[3]);
        }
        String filename = "normal_modes" + D + "D_"+nA+"_"+((int)(density*100))+"_cubic";
        if (args.length > 0) {
            filename = args[0];
        }
        
        System.out.println("Running "
                + (D == 1 ? "1D" : (D == 3 ? "FCC" : "2D hexagonal"))
                + " hard sphere simulation");
        System.out.println(nA + " atoms at density " + density);
        System.out.println(nSteps + " time units");
        System.out.println("output data to " + filename);

        // construct simulation
        SimDegreeFreedom sim = new SimDegreeFreedom(Space.getInstance(D), nA, density);
        
        // start simulation
        sim.activityIntegrate.setMaxSteps(nSteps/10);
        sim.getController().actionPerformed();
        System.out.println("equilibration finished");
        sim.getController().reset();
       
        sim.activityIntegrate.setMaxSteps(nSteps);
        sim.getController().actionPerformed();
        
        
        
        /* loops that
         *  -changes filename
         *  -changes accumulator histogram
         *  -changes histogram
         *  -calls actionPerformed
         *  
         */
        int accumulatorLength = sim.hists.length;
        for(int i = 0; i < accumulatorLength; i++){
            
            WriteHistograms wh;
            String outputName = new String("hist_" + i);
            wh = new WriteHistograms(outputName);
            wh.setHistogram(sim.hists[i].getHistograms());
            wh.actionPerformed();
        }
        
        
        
        
        
        System.out.println("Fini.");
    }

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
    MeterNormalModeCoordinate meternmc;
    WaveVectorFactory waveVectorFactory;
    MCMoveAtomCoupled mcmove;
    AccumulatorHistogram[] hists;

}