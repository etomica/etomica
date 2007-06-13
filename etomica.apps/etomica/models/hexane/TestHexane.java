/*
 * Created on May 24, 2005
 */
package etomica.models.hexane;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.BravaisLattice;
import etomica.lattice.crystal.Primitive;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.MCMoveMoleculeCoupled;
import etomica.phase.Phase;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
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

public class TestHexane extends Simulation {

    public ActivityIntegrate activityIntegrate;
    public IntegratorMC integrator;

    public Phase phase;

    public BoundaryDeformablePeriodic bdry;
    public BravaisLattice lattice;
    public CoordinateDefinition coordinateDefinition;
    public Primitive primitive;
    
//    public MCMoveVolume moveVolume;
//    public MCMoveCrankshaft crank; 
//    public MCMoveReptate snake;
    public MCMoveMolecule moveMolecule;
    public CBMCGrowSolidHexane growMolecule;
    public MCMoveRotateMolecule3D rot;
    public MCMoveMoleculeCoupled coupledMove;
    public MCMoveCombinedCbmcTranslation cctMove;

//    public PairIndexerMolecule pri;

    
    public TestHexane(Space space, int numMolecules, double dens) {
        //super(space, false, new PotentialMasterNbr(space, 12.0));
//        super(space, true, new PotentialMasterList(space, 12.0));
        super(space, false);
        PotentialMaster potentialMaster = new PotentialMaster(this);
        int chainLength = 6;
        int numAtoms = numMolecules * chainLength;
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
        defaults.makeLJDefaults();
        defaults.atomSize = 1.0;
        defaults.ignoreOverlap = false;

        SpeciesHexane species = new SpeciesHexane(this);
        getSpeciesManager().addSpecies(species);
        int[] nCells = new int[]{4,6,6};
        bdry = new BoundaryDeformableLattice(primitive, getRandom(), nCells);
        phase = new Phase(bdry);
        addPhase(phase);
        phase.getAgent(species).setNMolecules(numMolecules);
//        config.initializeCoordinates(phase);
        integrator = new IntegratorMC(potentialMaster, getRandom(),
                defaults.temperature);
        
        moveMolecule = new MCMoveMolecule(potentialMaster, getRandom(),
                0.1, 1, false);
        // 0.025 for translate, 0.042 for rotate for rho=0.3737735
        moveMolecule.setStepSize(0.024);        
        integrator.getMoveManager().addMCMove(moveMolecule);
        ((MCMoveStepTracker)moveMolecule.getTracker()).setNoisyAdjustment(true);
        
        // moveVolume = new MCMoveVolume(potentialMaster, phase.space(),
        // sim.getDefaults().pressure);
        // moveVolume.setPhase(phase);
        // integrator.getMoveManager().addMCMove(moveVolume);
        
        // crank = new MCMoveCrankshaft();

        // snake = new MCMoveReptate(this);
        // snake.setPhase(phase);
        // integrator.getMoveManager().addMCMove(snake);
        
        rot = new MCMoveRotateMolecule3D(potentialMaster, getRandom());
        rot.setPhase(phase);
        rot.setStepSize(0.042);
        integrator.getMoveManager().addMCMove(rot);
        ((MCMoveStepTracker)rot.getTracker()).setNoisyAdjustment(true);
        
        growMolecule = new CBMCGrowSolidHexane(potentialMaster,
                getRandom(), integrator, phase, species, 20);
        growMolecule.setPhase(phase);
        integrator.getMoveManager().addMCMove(growMolecule);

        coupledMove = new MCMoveMoleculeCoupled(potentialMaster, getRandom());
        integrator.getMoveManager().addMCMove(coupledMove);
        
        cctMove = new MCMoveCombinedCbmcTranslation(potentialMaster, growMolecule, getRandom());
        cctMove.setPhase(phase);
        integrator.getMoveManager().addMCMove(cctMove);
        
        // nan we're going to need some stuff in there to set the step sizes and
        // other stuff like that.

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
//        //		The atoms must be on the same molecule- CriterionMolecular
//        //		The atoms must be separated by 3 bonds, or 2 other atoms.
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
        coordinateDefinition = new CoordinateDefinitionHexane(phase, primitive, species);
        coordinateDefinition.initializeCoordinates(nCells);

        integrator.setPhase(phase);
        
        //nan this will need to be changed
//        pri = new PairIndexerMolecule(phase, new PrimitiveHexane(space));
    }

    public static void main(String[] args) {
        int numMolecules = 144; //144
        long nSteps = 10000;
        // Monson reports data for 0.373773507616 and 0.389566754417
        double density = 0.373773507616;

        boolean graphic = true;

        //parse arguments
        if(args.length > 1){
            nSteps = Long.parseLong(args[1]);
        }
        if(args.length > 2){
            numMolecules = Integer.parseInt(args[2]);
        }
        if(args.length > 3){
            density = Double.parseDouble(args[3]);
            if(density == 0.37) {density = 0.373773507616;}
            if(density == 0.40) {density = 0.389566754417;}
        }
        
        //spaces are now singletons; we can only have one instance, so we call
        // it with this method, not a "new" thing.
        TestHexane sim = new TestHexane(Space3D.getInstance(), numMolecules, density);

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
//            meterNormalMode.setPhase(sim.phase);

//            PhaseInflateDeformable pid = new PhaseInflateDeformable(sim.getSpace());
////            PhaseInflate pid = new PhaseInflate(sim.phase);
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
//            normalModeData.TE(1.0/(sim.phase.getSpeciesMaster().moleculeCount()*meterNormalMode.getCallCount()));
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