package etomica.models.hexane;

import java.util.ArrayList;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeGroup;
import etomica.atom.AtomTypeSphere;
import etomica.data.AccumulatorAverage;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.BravaisLattice;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.lattice.crystal.PrimitiveFcc;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.BoltzmannProcessor;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinitionLeaf;
import etomica.normalmode.MCMoveHarmonic;
import etomica.normalmode.MCMoveMoleculeCoupled;
import etomica.normalmode.MeterHarmonicEnergy;
import etomica.normalmode.MeterNormalMode;
import etomica.normalmode.NormalModes;
import etomica.normalmode.NormalModes1DHR;
import etomica.normalmode.NormalModesFromFile;
import etomica.normalmode.P2XOrder;
import etomica.normalmode.SimHarmonic;
import etomica.normalmode.WaveVectorFactory;
import etomica.box.Box;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential;
import etomica.potential.Potential2;
import etomica.potential.Potential2HardSpherical;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Pixel;

public class SimHarmonicHexane extends Simulation {

    public ActivityIntegrate activityIntegrate;
    public IntegratorMC integrator;

    public Box box;

    public BoundaryDeformablePeriodic bdry;
    public BravaisLattice lattice;
    public CoordinateDefinition coordinateDefinition;
    public Primitive primitive;
    
    public MCMoveMolecule moveMolecule;
    public CBMCGrowSolidHexane growMolecule;
    public MCMoveRotateMolecule3D rot;
    public MCMoveMoleculeCoupled coupledMove;
    public MCMoveCombinedCbmcTranslation cctMove;

//    public PairIndexerMolecule pri;

    
    public SimHarmonicHexane(Space space, double dens, int xCells, int yCells, int zCells) {
        //super(space, false, new PotentialMasterNbr(space, 12.0));
//        super(space, true, new PotentialMasterList(space, 12.0));
        super(space, false);
        PotentialMaster potentialMaster = new PotentialMaster(space);
        int chainLength = 6;
        //One molecule per cell
        int numAtoms = xCells * yCells * zCells * chainLength;
        primitive = new PrimitiveHexane(space);
        // close packed density is 0.4165783882178116
        // Monson reports data for 0.373773507616 and 0.389566754417
        primitive.scaleSize(Math.pow(0.4165783882178116/dens,1.0/3.0));
        lattice = new BravaisLattice(primitive);

        //This is the factor that multiples by the range of the potential in
        // order to define the area/volume in which neighbors are searched for.
        //This becomes the bond delta, which is the percentage the bond can
        // stretch, and I assume compress.
        double neighborRangeFac = 1.2;

        double bondFactor = 0.4;

        SpeciesHexane species = new SpeciesHexane(this);
        getSpeciesManager().addSpecies(species);
        int[] nCells = new int[]{xCells, yCells, zCells};
        bdry = new BoundaryDeformableLattice(primitive, getRandom(), nCells);
        box = new Box(bdry);
        addBox(box);
        box.getAgent(species).setNMolecules(xCells * yCells * zCells);
//        config.initializeCoordinates(box);
        integrator = new IntegratorMC(potentialMaster, getRandom(), 1.0);
        
        moveMolecule = new MCMoveMolecule(potentialMaster, getRandom(),
                0.1, 1, false);
        // 0.025 for translate, 0.042 for rotate for rho=0.3737735
        moveMolecule.setStepSize(0.024);        
        integrator.getMoveManager().addMCMove(moveMolecule);
        ((MCMoveStepTracker)moveMolecule.getTracker()).setNoisyAdjustment(true);
        
        // moveVolume = new MCMoveVolume(potentialMaster, box.space(),
        // sim.getDefaults().pressure);
        // moveVolume.setBox(box);
        // integrator.getMoveManager().addMCMove(moveVolume);
        
        // crank = new MCMoveCrankshaft();

        // snake = new MCMoveReptate(this);
        // snake.setBox(box);
        // integrator.getMoveManager().addMCMove(snake);
        
        rot = new MCMoveRotateMolecule3D(potentialMaster, getRandom());
        rot.setBox(box);
        rot.setStepSize(0.042);
        integrator.getMoveManager().addMCMove(rot);
        ((MCMoveStepTracker)rot.getTracker()).setNoisyAdjustment(true);
        
        growMolecule = new CBMCGrowSolidHexane(potentialMaster,
                getRandom(), integrator, box, species, 20);
        growMolecule.setBox(box);
        integrator.getMoveManager().addMCMove(growMolecule);

        coupledMove = new MCMoveMoleculeCoupled(potentialMaster, getRandom());
        integrator.getMoveManager().addMCMove(coupledMove);
        
        cctMove = new MCMoveCombinedCbmcTranslation(potentialMaster, growMolecule, getRandom());
        cctMove.setBox(box);
        integrator.getMoveManager().addMCMove(cctMove);
        
        // nan we're going to need some stuff in there to set the step sizes and
        // other stuff like that.

        integrator.setIsothermal(true);
        activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setMaxSteps(2000000);
        getController().addAction(activityIntegrate);
            
        //nan The box size we want is 5.72906360610622 by 11.21417818673970 by
        // 7.30591061708510
        //nan this is where the squared, unsquared box stuff comes in.
        //makes the density 0.41657 per Dr. Monson's comment in e-mail.
//        defaults.boxSize = 7.018;
//        defaults.boxSize = 100;

        //INTERMOLECULAR POTENTIAL STUFF

        //This potential is the intermolecular potential between atoms on
        // different molecules. We use the class "Potential" because we are
        // reusing the instance as we define each potential.
        Potential potential = new P2HardSphere(space);
        
        //here, we add the species to the PotentialMaster, using types.
        //The PotentialMaster generates a group potential and automatically
        // does a lot of the stuff which we have to do for the intramolecular
        // potential manually.
        AtomTypeSphere sphereType = (AtomTypeSphere) ((AtomTypeGroup) species
                .getMoleculeType()).getChildTypes()[0];

        //Add the Potential to the PotentialMaster
        potentialMaster.addPotential(potential, new AtomType[] { sphereType,
                sphereType });
        
        coupledMove.setPotential(potentialMaster.getPotential(new AtomType[] {
                species.getMoleculeType(), species.getMoleculeType() }  ));

        
        
//         //INTRAMOLECULAR POTENTIAL STUFF
//
//        //This PotentialGroup will hold all the intramolecular potentials.
//        //We give 1 as the argument because we are using 1 molecule to iterate
//        // on. The actual interactions between the atoms on the molecules will
//        // be calculated by a Potential2, but their summation is the molecule's
//        //effect on itself, which is a Potential1, or a Potential with nBody =
//        // 1.
//        PotentialGroup potentialChainIntra = potentialMaster.makePotentialGroup(1);
//
//            //BONDED INTERACTIONS
//
//        // This potential simulates the bonds between atoms in a molecule.
//        // XXX It will be superceded by a set of MC moves at some point in the
//        // future.
//        //This potential uses hard sphere interactions to model the bonded
//        // interactions of the atoms of the molecule.
//        //We make the bonding length 0.4 * sigma per Malanoski 1999.
//        potential = new P2HardSphere(space, defaults.atomSize * bondFactor, 
//                defaults.ignoreOverlap);
//        
//        //We will need an atom pair iterator (Api) that runs through the atoms
//        // on a single molecule.
//        //The atom pair iterator (Api) runs through the atoms on a single
//        // molecule.
//        //  It has an inner loop and an outer loop.
//        ApiIntragroup bonded = ApiBuilder.makeAdjacentPairIterator();
//        //We add the Potential and its Iterator to the PotentialGroup, in one
//        // fell swoop. Yay us!
//        potentialChainIntra.addPotential(potential, bonded);
//        
//            //NONBONDED INTERACTIONS
//        //This potential describes the basic hard sphere interactions between
//        // 2 atoms of a molecule.
//        
//        //Only the atoms next to each other interact, so we have two criteria:
//        //        The atoms must be on the same molecule- CriterionMolecular
//        //        The atoms must be separated by 3 bonds, or 2 other atoms.
//        ApiIntragroup nonbonded = ApiBuilder.makeNonAdjacentPairIterator(2);
//        potentialChainIntra.addPotential(potential, nonbonded);
//        
//        potentialMaster.addPotential(potentialChainIntra, new AtomType[] { species.getMoleculeType() } );

        // //INTRAMOLECULAR POTENTIAL STUFF
        //
        // //This PotentialGroup will hold all the intramolecular potentials.
        // //We give 1 as the argument because we are using 1 molecule to
        // iterate
        // // on. The actual interactions between the atoms on the molecules
        // will
        // // be calculated by a Potential2, but their summation is the
        // molecule's
        // //effect on itself, which is a Potential1, or a Potential with nBody
        // =
        // // 1.
        // PotentialGroup potentialChainIntra =
        // potentialMaster.makePotentialGroup(1);
        //
        // //BONDED INTERACTIONS
        //
        // // This potential simulates the bonds between atoms in a molecule.
        // // XXX It will be superceded by a set of MC moves at some point in
        // the
        // // future.
        // //This potential uses hard sphere interactions to model the bonded
        // // interactions of the atoms of the molecule.
        // //We make the bonding length 0.4 * sigma per Malanoski 1999.
        // potential = new P2HardSphere(space, defaults.atomSize * bondFactor,
        // defaults.ignoreOverlap);
        //        
        // //We will need an atom pair iterator (Api) that runs through the
        // atoms
        // // on a single molecule.
        // //The atom pair iterator (Api) runs through the atoms on a single
        // // molecule.
        // // It has an inner loop and an outer loop.
        // ApiIntragroup bonded = ApiBuilder.makeAdjacentPairIterator();
        // //We add the Potential and its Iterator to the PotentialGroup, in one
        // // fell swoop. Yay us!
        // potentialChainIntra.addPotential(potential, bonded);
        //        
        // //NONBONDED INTERACTIONS
        // //This potential describes the basic hard sphere interactions between
        // // 2 atoms of a molecule.
        //        
        // //Only the atoms next to each other interact, so we have two
        // criteria:
        // // The atoms must be on the same molecule- CriterionMolecular
        // // The atoms must be separated by 3 bonds, or 2 other atoms.
        // ApiIntragroup nonbonded = ApiBuilder.makeNonAdjacentPairIterator(2);
        // potentialChainIntra.addPotential(potential, nonbonded);
        //        
        // potentialMaster.addPotential(potentialChainIntra, new AtomType[] {
        // species.getMoleculeType() } );


        //Initialize the positions of the atoms.
        coordinateDefinition = new CoordinateDefinitionHexane(box, primitive, species);
        coordinateDefinition.initializeCoordinates(nCells);

        integrator.setBox(box);
        
        //nan this will need to be changed
//        pri = new PairIndexerMolecule(box, new PrimitiveHexane(space));
    }

    public static void main(String[] args) {
//        int numMolecules = 144; //144

        int xLng = 4;
        int yLng = 4;
        int zLng = 6;
        
        long nSteps = 10000;
        // Monson reports data for 0.373773507616 and 0.389566754417
        double density = 0.373773507616;

        boolean graphic = false;

        //parse arguments
        //filename is element 0
        if(args.length > 1){
            nSteps = Long.parseLong(args[1]);
        }
        if(args.length > 2){
            density = Double.parseDouble(args[2]);
            if(density == 0.37) {density = 0.373773507616;}
            if(density == 0.40) {density = 0.389566754417;}
        }
        if(args.length > 5){
            xLng = Integer.parseInt(args[3]);
            yLng = Integer.parseInt(args[4]);
            zLng = Integer.parseInt(args[5]);
        }
        
        //spaces are now singletons; we can only have one instance, so we call
        // it with this method, not a "new" thing.
        TestHexane sim = new TestHexane(Space3D.getInstance(), density, xLng, yLng, zLng);

        System.out.println("Happy Goodness!!");

        if (graphic) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim);
            simGraphic.makeAndDisplayFrame();
        } else {

            String filename = "normal_modes_hexane";

            PrimitiveHexane primitive = (PrimitiveHexane)sim.lattice.getPrimitive();
//            // primitive doesn't need scaling.  The boundary was designed to be commensurate with the primitive
//            WaveVectorFactorySimple waveVectorFactory = new WaveVectorFactorySimple(primitive);
//            // we need to set this up now even though we don't use it during equilibration so that
//            // the meter can grab the lattice points
//            MeterNormalMode meterNormalMode = new MeterNormalMode();
//            meterNormalMode.setWaveVectorFactory(waveVectorFactory);
//            meterNormalMode.setCoordinateDefinition(sim.coordinateDefinition);
//            meterNormalMode.setBox(sim.box);

//            BoxInflateDeformable pid = new BoxInflateDeformable(sim.getSpace());
////            BoxInflate pid = new BoxInflate(sim.box);
//            MeterPressureByVolumeChange meterPressure = new MeterPressureByVolumeChange(sim.getSpace(), pid);
//            meterPressure.setIntegrator(sim.integrator);
//            AccumulatorAverage pressureAccumulator = new AccumulatorAverage(sim);
//            DataPump pressureManager = new DataPump(meterPressure, pressureAccumulator);
//            pressureAccumulator.setBlockSize(50);
//            new IntervalActionAdapter(pressureManager, sim.integrator);

            sim.activityIntegrate.setMaxSteps(nSteps/10);
            sim.getController().actionPerformed();
            System.out.println("equilibration finished");

//            ((MCMoveStepTracker)sim.moveMolecule.getTracker()).setTunable(false);
//            ((MCMoveStepTracker)sim.rot.getTracker()).setTunable(false);
            
            sim.getController().reset();
            sim.activityIntegrate.setMaxSteps(nSteps);
            
//            IntervalActionAdapter adapter = new IntervalActionAdapter(meterNormalMode);
//            adapter.setActionInterval(100);
//            sim.integrator.addListener(adapter);

            sim.getController().actionPerformed();
            
//            DataGroup normalModeData = (DataGroup)meterNormalMode.getData();
//            normalModeData.TE(1.0/(sim.box.getSpeciesMaster().moleculeCount()*meterNormalMode.getCallCount()));
//            int normalDim = meterNormalMode.getCoordinateDefinition().getCoordinateDim();
//            
//            IVector[] waveVectors = waveVectorFactory.getWaveVectors();
//            double[] coefficients = waveVectorFactory.getCoefficients();
//            
//            try {
//                FileWriter fileWriterQ = new FileWriter(filename+".Q");
//                FileWriter fileWriterS = new FileWriter(filename+".S");
//                for (int i=0; i<waveVectors.length; i++) {
//                    fileWriterQ.write(Double.toString(coefficients[i]));
//                    for (int j=0; j<waveVectors[i].getD(); j++) {
//                        fileWriterQ.write(" "+waveVectors[i].x(j));
//                    }
//                    fileWriterQ.write("\n");
//                    DataDoubleArray dataS = (DataDoubleArray)normalModeData.getData(i);
//                    for (int k=0; k<normalDim; k++) {
//                        fileWriterS.write(Double.toString(dataS.getValue(k*normalDim)));
//                        for (int l=1; l<normalDim; l++) {
//                            fileWriterS.write(" "+dataS.getValue(k*normalDim+l));
//                        }
//                        fileWriterS.write("\n");
//                    }
//                }
//                fileWriterQ.close();
//                fileWriterS.close();
//            }
//            catch (IOException e) {
//                throw new RuntimeException("Oops, failed to write data "+e);
//            }
            
//            double avgPressure = ((DataDoubleArray)(((DataGroup)pressureAccumulator.getData()).getData(StatType.AVERAGE.index))).x;
//              avgPressure = ((DataDoubleArray)((DataGroup)pressureAccumulator.getData()).getData(AccumulatorAverage.StatType.AVERAGE.index)).x;
//              System.out.println("Avg Pres = "+ avgPressure);
        }

        System.out.println("Go look at the data!");
    }
}
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
//    private static final String APP_NAME = "Sim Harmonic";
//    private static final long serialVersionUID = 1L;
//    public IntegratorMC integrator;
//    public ActivityIntegrate activityIntegrate;
//    public Box box;
//    public Boundary boundary;
//    public Species species;
//    public NormalModes normalModes;
//    public int[] nCells;
//    public CoordinateDefinition coordinateDefinition;
//    public PrimitiveHexane primitive;
//    
//    private double density;
//   
//    
//    public SimHarmonicHexane(Space space, double dens, int xCells, int yCells, 
//            int zCells, String filename, double harmonicFudge) {
//        super(space, true);
//
//        int chainLength = 6;
//        int numAtoms = xCells * yCells * zCells * chainLength;
//        int D = space.D();
//        if( D != 3) {
//            throw new IllegalStateException("SimHarmonicHexane requires a 3D space");
//        }
//
//        species = new SpeciesHexane(this);
//        getSpeciesManager().addSpecies(species);
//
//        box = new Box(this);
//        addBox(box);
//        box.getAgent(species).setNMolecules(numAtoms);
//
//        integrator = new IntegratorMC(this, null);
//
//        activityIntegrate = new ActivityIntegrate(integrator);
//        getController().addAction(activityIntegrate);
//
//        MCMoveHarmonic move = new MCMoveHarmonic(getRandom());
//        integrator.getMoveManager().addMCMove(move);
//
//            primitive = new PrimitiveHexane(space);
//            double v = primitive.unitCell().getVolume();
//            // close packed density is 0.4165783882178116
//            // Monson reports data for 0.373773507616 and 0.389566754417
//            
//            primitive.scaleSize(Math.pow(v*dens,-1.0/3.0));
//            int n = (int)Math.round(Math.pow(numAtoms, 1.0/3.0));
//            nCells = new int[]{xCells, yCells, zCells};
//            boundary = new BoundaryDeformableLattice(primitive, getRandom(), nCells);
//
//            box.setBoundary(boundary);
//
//        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive);
//        coordinateDefinition.initializeCoordinates(nCells);
//        
//        if(D == 1) {
//            normalModes = new NormalModes1DHR();
//        } else {
//            normalModes = new NormalModesFromFile(filename, D);
//        }
//        normalModes.setHarmonicFudge(harmonicFudge);
//        
//        WaveVectorFactory waveVectorFactory = normalModes.getWaveVectorFactory();
//        waveVectorFactory.makeWaveVectors(box);
//        move.setOmegaSquared(normalModes.getOmegaSquared(box), waveVectorFactory.getCoefficients());
//        move.setEigenVectors(normalModes.getEigenvectors(box));
//        move.setWaveVectors(waveVectorFactory.getWaveVectors());
//        move.setWaveVectorCoefficients(waveVectorFactory.getCoefficients());
//        move.setCoordinateDefinition(coordinateDefinition);
//        
//        move.setBox(box);
//        
//        integrator.setBox(box);
//    }
//
//    /**
//     * @param args
//     */
//    public static void main(String[] args) {
//        
//        //set up simulation parameters
//        int D = 3;
//        int nA = 27;
//        
//        double density = 1.3;
//        double harmonicFudge = 1;
//        long steps = 800000;
//        if (D == 1) {
//            nA = 3;
//            density = 0.5;
//        }
//        boolean graphic = false;
//        String filename = "normal_modes3D_27_130";
//        if (args.length > 0) {
//            filename = args[0];
//        }
//        if (args.length > 1) {
//            density = Double.parseDouble(args[1]);
//        }
//        if (args.length > 2) {
//            steps = Long.parseLong(args[2]);
//        }
//        if (args.length > 3) {
//            nA = Integer.parseInt(args[3]);
//        }
//        if (args.length > 4) {
//            harmonicFudge = Double.parseDouble(args[4]);
//        }
//        System.out.println("Running "+(D==1 ? "1D" : (D==3 ? "FCC" : "2D hexagonal")) +" harmonic simulation, measuring hard sphere energy");
//        System.out.println(nA+" atoms at density "+density);
//        System.out.println("harmonic fudge: "+harmonicFudge);
//        System.out.println(steps+" MC steps");
//        
//        //construct simulation
//        SimHarmonic sim = new SimHarmonic(Space.getInstance(D), nA, density, filename, harmonicFudge);
//        
//        //add hard potentials for FEP calculations.  With de novo sampling potential is not otherwise used.
//        Potential2 p2 = new P2HardSphere(sim.getSpace(), 1.0, true);
//        if (D == 1) {
//            p2 = new P2XOrder(sim.getSpace(), (Potential2HardSpherical)p2);
//        }
//        PotentialMaster potentialMaster = (D == 1 ? new PotentialMasterList(sim) : new PotentialMaster(sim.getSpace()));
//        potentialMaster.addPotential(p2, new AtomType[]{sim.species.getMoleculeType(),sim.species.getMoleculeType()});
//
//        if (potentialMaster instanceof PotentialMasterList) {
//            double neighborRange;
//            if (D == 1) {
//                neighborRange = 1.01 / density;
//            }
//            else {
//                //FCC
//                double L = Math.pow(0.26*density, 1.0/3.0);
//                neighborRange = L / Math.sqrt(2.0);
//            }
//            ((PotentialMasterList)potentialMaster).setRange(neighborRange);
//            // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
//            ((PotentialMasterList)potentialMaster).getNeighborManager(sim.box).reset();
//        }
//
//        //meters for FEP calculations
//        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster);
//        meterPE.setBox(sim.box);
//        BoltzmannProcessor bp = new BoltzmannProcessor();
//        bp.setTemperature(1);
//        DataPump pump = new DataPump(meterPE,bp);
//        AccumulatorAverage avgBoltzmann = new AccumulatorAverage(1);
//        bp.setDataSink(avgBoltzmann);
//        avgBoltzmann.setPushInterval(5);
//        sim.integrator.addIntervalAction(pump);
//
////         MeterMomentumCOM meterCOM = new MeterMomentumCOM(sim.space);
////         MeterPositionCOM meterCOM = new MeterPositionCOM(sim.space);
////         DataSinkConsole console = new DataSinkConsole();
////         DataProcessorFunction filter = new DataProcessorFunction(new Function.Chop());
////         DataPump comPump = new DataPump(meterCOM,filter);
////         filter.setDataSink(console);
////         IntervalActionAdapter comAdapter = new IntervalActionAdapter(comPump);
////         sim.integrator.addListener(comAdapter);
////         meterCOM.setBox(sim.box);
//
//        //set up things for determining energy of harmonic system
//        //read and set up wave vectors
//
//        if(graphic){
//            //meter for harmonic system energy, sent to direct and to boltzmann average
//            MeterHarmonicEnergy harmonicEnergy = new MeterHarmonicEnergy(sim.coordinateDefinition, sim.normalModes);
//            harmonicEnergy.setBox(sim.box);
//            DataFork harmonicFork = new DataFork();
//            AccumulatorAverage harmonicAvg = new AccumulatorAverage(5);
//            DataPump pumpHarmonic = new DataPump(harmonicEnergy, harmonicFork);
//            harmonicFork.addDataSink(harmonicAvg);
//            sim.integrator.addIntervalAction(pumpHarmonic);
//
//            //histogram energy of individual modes
////            MeterHarmonicSingleEnergy harmonicSingleEnergy = new MeterHarmonicSingleEnergy(coordinateDefinitionLeaf, sim.normalModes);
////            harmonicSingleEnergy.setTemperature(1.0);
////            harmonicSingleEnergy.setBox(sim.box);
////    //        DataProcessorFunction harmonicLog = new DataProcessorFunction(new Function.Log());
////            AccumulatorAverage harmonicSingleAvg = new AccumulatorAverage(5);
////            DataHistogram harmonicSingleHistogram = new DataHistogram(new HistogramSimple.Factory(50, new DoubleRange(0, 1)));
////            pump = new DataPump(harmonicSingleEnergy, harmonicSingleHistogram);
////    //        harmonicLog.setDataSink(harmonicSingleHistogram);
////            harmonicSingleHistogram.setDataSink(harmonicSingleAvg);
////            iaa= new IntervalActionAdapter(pump);
////            iaa.setActionInterval(1);
////            sim.integrator.addListener(iaa);
//            
//            //set up measurement of S matrix, to check that configurations are generated as expected
//            MeterNormalMode meterNormalMode = new MeterNormalMode();
//            meterNormalMode.setCoordinateDefinition(sim.coordinateDefinition);
//            WaveVectorFactory waveVectorFactory = sim.normalModes.getWaveVectorFactory();
//            meterNormalMode.setWaveVectorFactory(waveVectorFactory);
//            meterNormalMode.setBox(sim.box);
//
//
//            //graphic simulation -- set up window
////            sim.getDefaults().pixelUnit = new Pixel(0.05);
//            SimulationGraphic simG = new SimulationGraphic(sim, APP_NAME);
//            ArrayList dataStreamPumps = simG.getController().getDataStreamPumps();
//            dataStreamPumps.add(pump);
//            dataStreamPumps.add(pumpHarmonic);
//            
//            DisplayTextBoxesCAE boxesPE = new DisplayTextBoxesCAE();
//            boxesPE.setAccumulator(avgBoltzmann);
//            boxesPE.setPrecision(6);
//            simG.add(boxesPE);
//
//            DisplayTextBoxesCAE harmonicBoxes = new DisplayTextBoxesCAE();
//            harmonicBoxes.setAccumulator(harmonicAvg);
//            simG.add(harmonicBoxes);
//
////            DisplayPlot harmonicPlot = new DisplayPlot();
////            harmonicPlot.setDoLegend(false);
////            harmonicSingleAvg.addDataSink(harmonicPlot.getDataSet().makeDataSink(), new StatType[]{StatType.AVERAGE});
////            simG.add(harmonicPlot);
//
//            simG.getDisplayBox(sim.box).setPixelUnit(new Pixel(10));
//            simG.makeAndDisplayFrame(APP_NAME);
//        } else {
//            //not graphic, so run simulation batch
//            //S data is written to file
//            sim.activityIntegrate.setMaxSteps(steps);
//
//            sim.getController().actionPerformed();
//            
//            DataGroup boltzmannData = (DataGroup)avgBoltzmann.getData();
//            double pNotOverlap = ((DataDouble)boltzmannData.getData(StatType.AVERAGE.index)).x;
//            double pError = ((DataDouble)boltzmannData.getData(StatType.ERROR.index)).x;
//            
//            System.out.println("avg HS Boltzmann factor "+pNotOverlap+" +/- "+pError);
//            
//            System.out.println("free energy contribution "+(-Math.log(pNotOverlap))+" +/- "+(pError/pNotOverlap));
//            System.out.println("free energy contribution per molecule "+(-Math.log(pNotOverlap)/nA)+" +/- "+(pError/pNotOverlap)/nA);
//        }
//
//    }
//
//
//}
