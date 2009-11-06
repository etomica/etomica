package etomica.models.oneDHardRods;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistogram;
import etomica.data.DataPump;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.listener.IntegratorListenerAction;
import etomica.math.SpecialFunctions;
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
import etomica.util.ReadParameters;

/**
 * MC simulation
 * 
 * 1D hard rods
 * No graphic display
 * Output: histogram files of probability that a mode is zero
 * Calculate free energy of solid
 * 
 * Treats modes as degrees of freedom; uses a MeterDifferentImage to calculate 
 * what happens when an extra mode is added.
 * 
 */

/*
 * Starts in notes 7/09
 */
public class SimDifferentImage1DHR extends Simulation {

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
    MeterDifferentImage meterdi;
    WaveVectorFactory waveVectorFactory;
    MCMoveAtomCoupled mcMoveAtom;
    MCMoveChangeMultipleWV mcMoveMode;
    AccumulatorHistogram[] hists;
    int harmonicWV;
    boolean[] skipThisMode;
    AccumulatorAverageFixed accumulatorDI;


    public SimDifferentImage1DHR(Space _space, int numAtoms, double density, 
            int blocksize, int[] changeable) {
        super(_space);
        
//        long seed = 3;
//        System.out.println("Seed explicitly set to " + seed);
//        IRandom rand = new RandomNumberGenerator(seed);
//        this.setRandom(rand);
        
        PotentialMasterList potentialMaster = new PotentialMasterList(this, space);

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);
        
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
        
        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(nCells);
        int coordinateDim = coordinateDefinition.getCoordinateDim();

        double neighborRange = 1.01/density;
        potentialMaster.setRange(neighborRange);
        //find neighbors now.  Don't hook up NeighborListManager since the
        //  neighbors won't change
        potentialMaster.getNeighborManager(box).reset();
        
        integrator = new IntegratorMC(this, potentialMaster);
        integrator.setBox(box);
        
        nm = new NormalModes1DHR(box.getBoundary(), numAtoms);
        nm.setHarmonicFudge(1.0);
        nm.setTemperature(1.0);
        nm.getOmegaSquared();
        waveVectorFactory = nm.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);
        
        //Set up skip-these-modes code
        double[] wvc= nm.getWaveVectorFactory().getCoefficients();
        double[][] omega = nm.getOmegaSquared();
        int jump = coordinateDim * nm.getWaveVectorFactory().getWaveVectors().length;
        skipThisMode = new boolean[2*jump];
        for(int i = 0; i < 2*jump; i++){
            skipThisMode[i] = false;
        }
        for(int wvCount = 0; wvCount < wvc.length; wvCount++){
            //Sets up the imaginary modes that should be skipped.
            if(wvc[wvCount] == 0.5) {
                for(int j = 0; j < coordinateDim; j++){
                    skipThisMode[j + coordinateDim*wvCount + jump] = true;
                }
            }
            //Sets up the modes that are center of mass motion to skip
            for(int j = 0; j < omega[wvCount].length; j++){
                if(Double.isInfinite(omega[wvCount][j])){
                    skipThisMode[j + coordinateDim*wvCount] = true;
                    skipThisMode[j + coordinateDim*wvCount + jump] = true;

                }
            }
        }
        
        mcMoveAtom = new MCMoveAtomCoupled(potentialMaster, random, space);
        mcMoveAtom.setPotential(potential);
        mcMoveAtom.setBox(box);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        mcMoveAtom.setStepSizeMin(0.001);
        mcMoveAtom.setStepSize(0.01);
        
        mcMoveMode = new MCMoveChangeMultipleWV(potentialMaster, random);
        mcMoveMode.setBox(box);
        integrator.getMoveManager().addMCMove(mcMoveMode);
        mcMoveMode.setCoordinateDefinition(coordinateDefinition);
        mcMoveMode.setEigenVectors(nm.getEigenvectors());
        mcMoveMode.setOmegaSquared(nm.getOmegaSquared());
        mcMoveMode.setWaveVectorCoefficients(nm.getWaveVectorFactory().getCoefficients());
        mcMoveMode.setWaveVectors(nm.getWaveVectorFactory().getWaveVectors());
        mcMoveMode.setChangeableWVs(changeable);
        
        meterdi = new MeterDifferentImage("MeterDI", potentialMaster, numAtoms, density, (Simulation)this, 
                primitive, basis, coordinateDefinition, nm);
        
        accumulatorDI = new AccumulatorAverageFixed(blocksize);
        DataPump pumpFromMeter = new DataPump(meterdi, accumulatorDI);
        
        IntegratorListenerAction pumpListener = new IntegratorListenerAction(pumpFromMeter);
        pumpListener.setInterval(blocksize);
        integrator.getEventManager().addListener(pumpListener);
        
        activityIntegrate = new ActivityIntegrate(integrator, 0, true);
        getController().addAction(activityIntegrate);
        
        
//        IAtomList leaflist = box.getLeafList();
//        double[] locations = new double[numAtoms];
//        System.out.println("starting positions:");
//        for(int i = 0; i < numAtoms; i++){
//            //one d is assumed here.
//            locations[i] = ( ((Atom)leaflist.getAtom(i)).getPosition().x(0) );
//        }
//        
//        for(int i = 0; i < numAtoms; i++){
//            System.out.println(i + "  " + locations[i]);
//        }
    }

    public Primitive getPrimitive() {
        return primitive;
    }

    public BasisMonatomic getBasis() {
        return basis;
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
        int[] changeableWV = params.changeableWV;
        long nSteps = params.numSteps;
        int bs = params.blockSize;
        String outputfn = params.outputname;
        
        System.out.println("Running "
                + (D == 1 ? "1D" : (D == 3 ? "FCC" : "2D hexagonal"))
                + " hard sphere simulation");
        System.out.println(nA + " atoms at density " + density);
        System.out.println(nSteps + " steps, " + bs + " blocksize");
        System.out.println("input data from " + inputFilename);
        System.out.println("output data to " + filename);

        // construct simulation
        SimDifferentImage1DHR sim = new SimDifferentImage1DHR(Space.getInstance(D),
                nA, density, bs, changeableWV);
        
        // start simulation
        sim.activityIntegrate.setMaxSteps(nSteps/10);
        sim.getController().actionPerformed();
        System.out.println("equilibration finished");
        sim.getController().reset();
       
        sim.activityIntegrate.setMaxSteps(nSteps);
        sim.getController().actionPerformed();
        
      //After processing...
        DataGroup group = (DataGroup)sim.accumulatorDI.getData();
        double results = ((DataDouble)group.getData(StatType.AVERAGE.index)).x;
        System.out.println("results: " + results);
        
//        IAtomList leaflist = sim.box.getLeafList();
//        double[] locations = new double[nA];
//        System.out.println("final:");
//        for(int i = 0; i < nA; i++){
//            //one d is assumed here.
//            locations[i] = ( ((Atom)leaflist.getAtom(i)).getPosition().x(0) );
//        }
//        
//        for(int i = 0; i < 32; i++){
//            System.out.println(i + "  " + locations[i]);
//        }
        
        if(D==1) {
            double AHR = -(nA-1)*Math.log(nA/density-nA)
                + SpecialFunctions.lnFactorial(nA) ;
            System.out.println("Hard-rod free energy: "+AHR);
        }
        
        System.out.println("Fini.");
    }
    
    public static class SimParam extends ParameterBase {
        public int numAtoms = 10;
        public double density = 0.70;
        public int D = 1;
        public double harmonicFudge = 1.0;
        public String filename = "HR1D_";
        public String inputfilename = "input";
        public String outputname = "hists";
        public double temperature = 1.0;
        public int[] changeableWV = {0, 1, 2, 3, 4};
        public int nBins = 200;
        
        public int blockSize = 1000;
        public long numSteps = 1000000;
    }

}