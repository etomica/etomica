/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.hexane;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.*;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.lattice.BravaisLattice;
import etomica.lattice.crystal.Primitive;
import etomica.listener.IntegratorListenerAction;
import etomica.normalmode.*;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.dimensions.Null;

import java.util.ArrayList;
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

	private static final String APP_NAME = "Test Hexane Harmonic";

    public TestHexaneHarmonic(Space _space, int numMolecules) {
        //super(space, false, new PotentialMasterNbr(space, 12.0));
//        super(space, true, new PotentialMasterList(space, 12.0));
        super(_space);
        PotentialMaster potentialMaster = new PotentialMaster();
        int chainLength = 6;
        int numAtoms = numMolecules * chainLength;
        primitive = new PrimitiveHexane(space);
        // close packed density is 0.4165783882178116
        // Monson reports data for 0.373773507616 and 0.389566754417
        primitive.scaleSize(Math.pow(0.4165783882178116/0.373773507616,1.0/3.0));
        lattice = new BravaisLattice(primitive);

        //This is the factor that multiples by the range of the potential in
        // order to define the area/volume in which neighbors are searched for.
        //This becomes the bond delta, which is the percentage the bond can
        // stretch, and I assume compress.
        double neighborRangeFac = 1.2;

        double bondFactor = 0.4;

        SpeciesHexane species = new SpeciesHexane(space);
        addSpecies(species);
        int[] nCells = new int[]{4,6,6};
        bdry =  new BoundaryDeformableLattice(primitive, nCells);
        box = new Box(bdry, space);
        addBox(box);
        box.setNMolecules(species, numMolecules);
//        config.initializeCoordinates(box);

        integrator = new IntegratorMC(potentialMaster, getRandom(), 1.0);
        moveMolecule = new MCMoveMolecule(potentialMaster, getRandom(), space, 1.0, 15.0);
//        moveVolume = new MCMoveVolume(potentialMaster, box.space(), sim.getDefaults().pressure);
//        moveVolume.setBox(box);
//        crank = new MCMoveCrankshaft();
        
//         snake = new MCMoveReptate(this);
//         snake.setBox(box);
         
         rot = new MCMoveRotateMolecule3D(potentialMaster, getRandom(), space);
         rot.setBox(box);
         
         // 0.025 for translate, 0.042 for rotate for rho=0.3737735
         moveMolecule.setStepSize(0.024);
         rot.setStepSize(0.042);

        //nan we're going to need some stuff in there to set the step sizes and other stuff like that.
        
        integrator.getMoveManager().addMCMove(moveMolecule);
//        integrator.getMoveManager().addMCMove(snake);
        integrator.getMoveManager().addMCMove(rot); 
//        integrator.getMoveManager().addMCMove(moveVolume);
        
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


        getController().addAction(activityIntegrate);

        //INTERMOLECULAR POTENTIAL STUFF

        //This potential is the intermolecular potential between atoms on
        // different molecules. We use the class "Potential" because we are
        // reusing the instance as we define each potential.
        Potential potential = new P2HardSphere(space);
        
        //here, we add the species to the PotentialMaster, using types.
        //The PotentialMaster generates a group potential and automatically
        // does a lot of the stuff which we have to do for the intramolecular
        // potential manually.
        AtomType sphereType = species.getLeafType();

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

        //Initialize the positions of the atoms.
        coordinateDefinition = new CoordinateDefinitionHexane(this, box, primitive, species, space);
        coordinateDefinition.initializeCoordinates(nCells);

        integrator.setBox(box);
        
        //nan this will need to be changed
//        pri = new PairIndexerMolecule(box, new PrimitiveHexane(space));
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
        
        NormalModes normalModes = new NormalModesFromFile(filename, 3);
        
        MeterHarmonicEnergy harmonicEnergy = new MeterHarmonicEnergy(sim.coordinateDefinition, normalModes);
        DataFork harmonicFork = new DataFork();
        AccumulatorAverageCollapsing harmonicAvg = new AccumulatorAverageCollapsing();
        DataPump pump = new DataPump(harmonicEnergy, harmonicFork);
        harmonicFork.addDataSink(harmonicAvg);
        IntegratorListenerAction pumpListener = new IntegratorListenerAction(pump);
        pumpListener.setInterval(100);
        sim.integrator.getEventManager().addListener(pumpListener);
        BoltzmannProcessor boltz = new BoltzmannProcessor();
        boltz.setTemperature(1.0);
        harmonicFork.addDataSink(boltz);
        AccumulatorAverageCollapsing harmonicBoltzAvg = new AccumulatorAverageCollapsing();
        boltz.setDataSink(harmonicBoltzAvg);
        DataProcessorFoo fooer = new DataProcessorFoo();
        harmonicBoltzAvg.addDataSink(fooer, new StatType[]{harmonicBoltzAvg.AVERAGE});
        
        MeterHarmonicSingleEnergy harmonicSingleEnergy = new MeterHarmonicSingleEnergy(sim.coordinateDefinition, normalModes);
        harmonicSingleEnergy.setBox(sim.box);
//        DataProcessorFunction harmonicLog = new DataProcessorFunction(new Function.Log());
        boltz = new BoltzmannProcessor();
        boltz.setTemperature(1);
        DataPump pumpSingle = new DataPump(harmonicSingleEnergy, boltz);
        //DataHistogram harmonicSingleHistogram = new DataHistogram(new HistogramSimple(50, new DoubleRange(0, 1)));
        AccumulatorAverageCollapsing harmonicSingleAvg = new AccumulatorAverageCollapsing();
        boltz.setDataSink(harmonicSingleAvg);
//        harmonicLog.setDataSink(harmonicSingleHistogram);
//        harmonicSingleHistogram.setDataSink(harmonicSingleAvg);
        //harmonicSingleAvg.addDataSink(harmonicSingleHistogram, new StatType[]{harmonicSingleAvg.AVERAGE});
        IntegratorListenerAction singlePumpListener = new IntegratorListenerAction(pumpSingle);
        singlePumpListener.setInterval(100);
        sim.integrator.getEventManager().addListener(singlePumpListener);
        DataProcessorFoo fooerSingle = new DataProcessorFoo();
        harmonicSingleAvg.addDataSink(fooerSingle, new StatType[]{harmonicSingleAvg.AVERAGE});

        if (graphic) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME, sim.space, sim.getController());
            ArrayList dataStreamPumps = simGraphic.getController().getDataStreamPumps();
            dataStreamPumps.add(pump);
            dataStreamPumps.add(pumpSingle);
            
            DisplayTextBoxesCAE harmonicBoxes = new DisplayTextBoxesCAE();
            harmonicBoxes.setAccumulator(harmonicAvg);
            simGraphic.add(harmonicBoxes);

//            DataHistogram harmonicSingleHistogram = new DataHistogram(new HistogramSimple.Factory(20, new DoubleRange(0, 1000)));
//            harmonicSingleAvg.addDataSink(harmonicSingleHistogram, new StatType[]{StatType.AVERAGE});
            DisplayPlot harmonicPlot = new DisplayPlot();
            harmonicPlot.setDoLegend(false);
            //harmonicSingleHistogram.setDataSink(harmonicPlot.getDataSet().makeDataSink());
            simGraphic.add(harmonicPlot);
            
            DisplayTextBox diffSingleA = new DisplayTextBox();
            diffSingleA.setLabel("deltaA, independent approx");
            fooerSingle.setDataSink(diffSingleA);
            simGraphic.add(diffSingleA);

            DisplayTextBox diffA = new DisplayTextBox();
            diffA.setLabel("deltaA");
            fooer.setDataSink(diffA);
            simGraphic.add(diffA);
            
            simGraphic.makeAndDisplayFrame(APP_NAME);
        }
        else {
            long nSteps = 10000;

            sim.activityIntegrate.setMaxSteps(nSteps);
            
            sim.getController().actionPerformed();

            double avgHarmonicEnergy = ((DataDouble)((DataGroup)harmonicAvg.getData()).getData(harmonicAvg.AVERAGE.index)).x;
            double errorHarmonicEnergy = ((DataDouble)((DataGroup)harmonicAvg.getData()).getData(harmonicAvg.ERROR.index)).x;
            System.out.println("avg harmonic energy: "+avgHarmonicEnergy+" +/- "+errorHarmonicEnergy);
        }

    }
    
    /**
     * DataProcessor that sums up the logs of all incoming values
     */
    public static class DataProcessorFoo extends DataProcessor {

        public DataPipe getDataCaster(IDataInfo incomingDataInfo) {
            return null;
        }
        
        public IDataInfo processDataInfo(IDataInfo incomingDataInfo) {
            dataInfo = new DataInfoDouble("free energy difference", Null.DIMENSION);
            data = new DataDouble();
            return dataInfo;
        }
            
        
        public IData processData(IData incomingData) {
            data.x = 0;
            int nData = incomingData.getLength();
            for (int i=0; i<nData; i++) {
                data.x += Math.log(incomingData.getValue(i));
            }
            return data;
        }
        
        private static final long serialVersionUID = 1L;
        private DataDouble data;
    }


    public ActivityIntegrate activityIntegrate;
    public IntegratorMC integrator;

    public Box box;

    public BoundaryDeformablePeriodic bdry;
 
    public MCMoveMolecule moveMolecule;
    public BravaisLattice lattice;
    public Primitive primitive;
    public CoordinateDefinition coordinateDefinition;
    
//    public MCMoveVolume moveVolume;
//    public MCMoveCrankshaft crank; 
//    public MCMoveReptate snake;
    
    public MCMoveRotateMolecule3D rot;
    
//    public PairIndexerMolecule pri;

}
