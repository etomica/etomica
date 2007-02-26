/*
 * Created on May 24, 2005
 */
package etomica.models.hexane;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.DataFork;
import etomica.data.DataHistogram;
import etomica.data.DataPump;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxesCAE;
import etomica.graphics.DisplayPhaseCanvas3DOpenGL;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntervalActionAdapter;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.lattice.BravaisLattice;
import etomica.normalmode.ArrayReader1D;
import etomica.normalmode.ArrayReader2D;
import etomica.normalmode.BoltzmannProcessor;
import etomica.normalmode.MeterHarmonicEnergy;
import etomica.normalmode.MeterHarmonicSingleEnergy;
import etomica.normalmode.MeterNormalMode;
import etomica.normalmode.NormalCoordLeaf;
import etomica.normalmode.WaveVectorFactorySimple;
import etomica.normalmode.TestFccHarmonic.DataProcessorFoo;
import etomica.phase.Phase;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.util.DoubleRange;
import etomica.util.HistogramSimple;
/**
 * @author nancycribbin
 *  
 */

/*
 * We use a PotentialMaster, rather than a PotentialMasterNbr, so that we do not
 * need to deal with cells, which BoundaryDeformablePeriodic cannot deal with at
 * this time.
 * 
 * @author nancycribbin
 *  
 */

public class TestHexaneHarmonic extends Simulation {

    public TestHexaneHarmonic(Space space, int numMolecules) {
        //super(space, false, new PotentialMasterNbr(space, 12.0));
//        super(space, true, new PotentialMasterList(space, 12.0));
        super(space, false, new PotentialMaster(space));
        int chainLength = 6;
        int numAtoms = numMolecules * chainLength;
        PrimitiveHexane primitive = new PrimitiveHexane(space);
        // close packed density is 0.4165783882178116
        // Monson reports data for 0.373773507616 and 0.389566754417
        primitive.scaleSize(Math.pow(0.4165783882178116/0.373773507616,1.0/3.0));
        lattice = new BravaisLattice(primitive);
        ConfigurationLattice config = new ConfigurationLattice(lattice);

        //This is the factor that multiples by the range of the potential in
        // order to define the area/volume in which neighbors are searched for.
        //This becomes the bond delta, which is the percentage the bond can
        // stretch, and I assume compress.
        double neighborRangeFac = 1.2;

        double bondFactor = 0.4;
        defaults.makeLJDefaults();
        defaults.atomSize = 1.0;
        defaults.ignoreOverlap = false;

        SpeciesHexane species = new SpeciesHexane(this);
        phase = new Phase(this);
        phase.getAgent(species).setNMolecules(numMolecules);
//        config.initializeCoordinates(phase);

        integrator = new IntegratorMC(potentialMaster, getRandom(), defaults.temperature);
        moveMolecule = new MCMoveMolecule(potentialMaster, getRandom(), defaults.atomSize, defaults.boxSize/2, false);
//        moveVolume = new MCMoveVolume(potentialMaster, phase.space(), sim.getDefaults().pressure);
//        moveVolume.setPhase(phase);
//        crank = new MCMoveCrankshaft();
        
//         snake = new MCMoveReptate(this);
//         snake.setPhase(phase);
         
         rot = new MCMoveRotateMolecule3D(getPotentialMaster(), getRandom());
         rot.setPhase(phase);
         
         // 0.025 for translate, 0.042 for rotate for rho=0.3737735
         moveMolecule.setStepSize(0.024);
         rot.setStepSize(0.042);

        //nan we're going to need some stuff in there to set the step sizes and other stuff like that.
        
        integrator.getMoveManager().addMCMove(moveMolecule);
//        integrator.getMoveManager().addMCMove(snake);
        integrator.getMoveManager().addMCMove(rot); 
//        integrator.getMoveManager().addMCMove(moveVolume);
        
        integrator.setIsothermal(true);
        activityIntegrate = new ActivityIntegrate(this, integrator);
        activityIntegrate.setMaxSteps(2000000);
        getController().addAction(activityIntegrate);
            
        //nan The box size we want is 5.72906360610622 by 11.21417818673970 by
        // 7.30591061708510
        //nan this is where the squared, unsquared box stuff comes in.
        //makes the density 0.41657 per Dr. Monson's comment in e-mail.
//        defaults.boxSize = 7.018;
//        defaults.boxSize = 100;


        getController().addAction(activityIntegrate);

        //INTERMOLECULAR POTENTIAL STUFF

        //This potential is the intermolecular potential between atoms on
        // different molecules. We use the class "Potential" because we are
        // reusing the instance as we define each potential.
        Potential potential = new P2HardSphere(space, defaults.atomSize, 
                defaults.ignoreOverlap);
        
        //here, we add the species to the PotentialMaster, using types.
        //The PotentialMaster generates a group potential and automatically
        // does a lot of the stuff which we have to do for the intramolecular
        // potential manually.
        AtomTypeSphere sphereType = (AtomTypeSphere) ((AtomFactoryHomo) species
                .moleculeFactory()).getChildFactory().getType();

        //Add the Potential to the PotentialMaster
        potentialMaster.addPotential(potential, new AtomType[] { sphereType,
                sphereType });
        
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
//        //		The atoms must be on the same molecule- CriterionMolecular
//        //		The atoms must be separated by 3 bonds, or 2 other atoms.
//        ApiIntragroup nonbonded = ApiBuilder.makeNonAdjacentPairIterator(2);
//        potentialChainIntra.addPotential(potential, nonbonded);
//        
//        potentialMaster.addPotential(potentialChainIntra, new AtomType[] { species.getMoleculeType() } );

        bdry =  new BoundaryDeformableLattice(primitive, new int[]{4,6,6});
        phase.setBoundary(bdry);

        //Initialize the positions of the atoms.
        config.initializeCoordinates(phase);

        integrator.setPhase(phase);
        
        //nan this will need to be changed
//        pri = new PairIndexerMolecule(phase, new PrimitiveHexane(space));
    }

    public static void main(String[] args) {
        int numMolecules = 144; //144
        boolean graphic = true;

        //spaces are now singletons; we can only have one instance, so we call
        // it with this method, not a "new" thing.
        TestHexaneHarmonic sim = new TestHexaneHarmonic(Space3D.getInstance(), numMolecules);

        String filename = "normal_modes_hexane";
        if (args.length > 0) {
            filename = args[0];
        }
        
        double[][] omegaSquared = ArrayReader1D.getFromFile(filename+".val");
        for (int i=0; i<omegaSquared.length; i++) {
            for (int j=0; j<omegaSquared[i].length; j++) {
                // omega is sqrt(kT)/sqrt(eigenvalue)
                omegaSquared[i][j] = 1/omegaSquared[i][j];
            }
        }
        double[][] waveVectorsAndCoefficients = ArrayReader1D.getFromFile(filename+".Q");
        IVector[] waveVectors = new IVector[waveVectorsAndCoefficients.length];
        double[] coefficients = new double[waveVectors.length];
        for (int i=0; i<waveVectors.length; i++) {
            coefficients[i] = waveVectorsAndCoefficients[i][0];
            waveVectors[i] = new Vector3D(waveVectorsAndCoefficients[i][1],
                    waveVectorsAndCoefficients[i][2],
                    waveVectorsAndCoefficients[i][3]);
        }
        double[][][] eigenvectors = ArrayReader2D.getFromFile(filename+".vec");

        NormalCoordHexane normalCoordHexane = new NormalCoordHexane();
        
        MeterHarmonicEnergy harmonicEnergy = new MeterHarmonicEnergy();
        harmonicEnergy.setEigenvectors(eigenvectors);
        harmonicEnergy.setOmegaSquared(omegaSquared);
        harmonicEnergy.setWaveVectors(waveVectors, coefficients);
        harmonicEnergy.setNormalCoordWrapper(normalCoordHexane);
        harmonicEnergy.setPhase(sim.phase);
        DataFork harmonicFork = new DataFork();
        AccumulatorAverage harmonicAvg = new AccumulatorAverage(sim);
        DataPump pump = new DataPump(harmonicEnergy, harmonicFork);
        harmonicFork.addDataSink(harmonicAvg);
        IntervalActionAdapter adapter = new IntervalActionAdapter(pump);
        adapter.setActionInterval(100);
        sim.integrator.addListener(adapter);
        BoltzmannProcessor boltz = new BoltzmannProcessor();
        boltz.setTemperature(1.0);
        harmonicFork.addDataSink(boltz);
        AccumulatorAverage harmonicBoltzAvg = new AccumulatorAverage(50);
        boltz.setDataSink(harmonicBoltzAvg);
        DataProcessorFoo fooer = new DataProcessorFoo();
        harmonicBoltzAvg.addDataSink(fooer, new StatType[]{StatType.AVERAGE});
        sim.register(harmonicEnergy, pump);
        
        MeterHarmonicSingleEnergy harmonicSingleEnergy = new MeterHarmonicSingleEnergy();
        harmonicSingleEnergy.setEigenvectors(eigenvectors);
        harmonicSingleEnergy.setOmegaSquared(omegaSquared);
        harmonicSingleEnergy.setWaveVectors(waveVectors, coefficients);
        harmonicSingleEnergy.setNormalCoordMapper(normalCoordHexane);
        harmonicSingleEnergy.setPhase(sim.phase);
        harmonicSingleEnergy.setTemperature(1.0);
//        DataProcessorFunction harmonicLog = new DataProcessorFunction(new Function.Log());
        DataHistogram harmonicSingleHistogram = new DataHistogram(new HistogramSimple.Factory(50, new DoubleRange(0, 1)));
        AccumulatorAverage harmonicSingleAvg = new AccumulatorAverage(sim);
        pump = new DataPump(harmonicSingleEnergy, harmonicSingleAvg);
//        harmonicLog.setDataSink(harmonicSingleHistogram);
//        harmonicSingleHistogram.setDataSink(harmonicSingleAvg);
        harmonicSingleAvg.addDataSink(harmonicSingleHistogram, new StatType[]{StatType.AVERAGE});
        adapter = new IntervalActionAdapter(pump);
        adapter.setActionInterval(100);
        sim.integrator.addListener(adapter);
        DataProcessorFoo fooerSingle = new DataProcessorFoo();
        harmonicSingleAvg.addDataSink(fooerSingle, new StatType[]{StatType.AVERAGE});
        sim.register(harmonicSingleEnergy, pump);

        if (graphic) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim);
            // don't spend all of our time repainting
//            ((DisplayPhaseCanvas3DOpenGL)simGraphic.getDisplayPhase(sim.phase).canvas).setAnimateFps(1);
            DisplayBoxesCAE harmonicBoxes = new DisplayBoxesCAE();
            harmonicBoxes.setAccumulator(harmonicAvg);
            simGraphic.add(harmonicBoxes);

//            DataHistogram harmonicSingleHistogram = new DataHistogram(new HistogramSimple.Factory(20, new DoubleRange(0, 1000)));
//            harmonicSingleAvg.addDataSink(harmonicSingleHistogram, new StatType[]{StatType.AVERAGE});
            DisplayPlot harmonicPlot = new DisplayPlot();
            harmonicPlot.setDoLegend(false);
            harmonicSingleHistogram.setDataSink(harmonicPlot.getDataSet().makeDataSink());
            simGraphic.add(harmonicPlot);
            
            DisplayBox diffSingleA = new DisplayBox();
            diffSingleA.setLabel("deltaA, independent approx");
            fooerSingle.setDataSink(diffSingleA);
            simGraphic.add(diffSingleA);

            DisplayBox diffA = new DisplayBox();
            diffA.setLabel("deltaA");
            fooer.setDataSink(diffA);
            simGraphic.add(diffA);
            
            simGraphic.makeAndDisplayFrame();
        }
        else {
            long nSteps = 10000;

            sim.activityIntegrate.setMaxSteps(nSteps);
            
            sim.getController().actionPerformed();

            double avgHarmonicEnergy = ((DataDouble)((DataGroup)harmonicAvg.getData()).getData(AccumulatorAverage.StatType.AVERAGE.index)).x;
            double errorHarmonicEnergy = ((DataDouble)((DataGroup)harmonicAvg.getData()).getData(AccumulatorAverage.StatType.ERROR.index)).x;
            System.out.println("avg harmonic energy: "+avgHarmonicEnergy+" +/- "+errorHarmonicEnergy);
        }

    }

    public ActivityIntegrate activityIntegrate;
    public IntegratorMC integrator;

    public Phase phase;

    public BoundaryDeformablePeriodic bdry;
 
    public MCMoveMolecule moveMolecule;
    public BravaisLattice lattice;
    
//    public MCMoveVolume moveVolume;
//    public MCMoveCrankshaft crank; 
//    public MCMoveReptate snake;
    
    public MCMoveRotateMolecule3D rot;
    
//    public PairIndexerMolecule pri;

}