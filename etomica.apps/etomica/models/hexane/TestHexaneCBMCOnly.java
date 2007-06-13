/*
 * Created on May 24, 2005
 */
package etomica.models.hexane;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.config.ConfigurationLattice;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.lattice.BravaisLattice;
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

public class TestHexaneCBMCOnly extends Simulation {

	private static final String APP_NAME = "Test Hexane CBMC Only";

    public TestHexaneCBMCOnly(Space space, int numMolecules) {
        // super(space, false, new PotentialMasterNbr(space, 12.0));
        // super(space, true, new PotentialMasterList(space, 12.0));
        super(space, false);
        PotentialMaster potentialMaster = new PotentialMaster(space);
        int chainLength = 6;
        int numAtoms = numMolecules * chainLength;
        PrimitiveHexane primitive = new PrimitiveHexane(space);
        // close packed density is 0.4165783882178116
        // Monson reports data for 0.373773507616 and 0.389566754417
        primitive.scaleSize(Math.pow(0.4165783882178116 / 0.373773507616,
                1.0 / 3.0));
//        primitive.scaleSize(Math.pow(0.4165783882178116 / 0.01,
//                1.0 / 3.0));
        lattice = new BravaisLattice(primitive);
        ConfigurationLattice config = new ConfigurationLattice(lattice);

        // This is the factor that multiples by the range of the potential in
        // order to define the area/volume in which neighbors are searched for.
        // This becomes the bond delta, which is the percentage the bond can
        // stretch, and I assume compress.
        double neighborRangeFac = 1.2;

        double bondFactor = 0.4;
        defaults.makeLJDefaults();
        defaults.atomSize = 1.0;
        defaults.ignoreOverlap = false;

        SpeciesHexane species = new SpeciesHexane(this);
        getSpeciesManager().addSpecies(species);
        bdry = new BoundaryDeformableLattice(primitive, getRandom(), new int[] {
            4, 6, 6 });
        phase = new Phase(bdry);
        addPhase(phase);
        phase.getAgent(species).setNMolecules(numMolecules);
        // config.initializeCoordinates(phase);

        integrator = new IntegratorMC(potentialMaster, getRandom(),
                defaults.temperature);

        growMolecule = new CBMCGrowSolidHexane(potentialMaster,
                getRandom(), integrator, phase, species, 20);
        growMolecule.setPhase(phase);
        integrator.getMoveManager().addMCMove(growMolecule);

        // nan we're going to need some stuff in there to set the step sizes and
        // other stuff like that.

        integrator.setIsothermal(true);
        activityIntegrate = new ActivityIntegrate(this, integrator);
        activityIntegrate.setMaxSteps(2000000);
        getController().addAction(activityIntegrate);

        // nan The box size we want is 5.72906360610622 by 11.21417818673970 by
        // 7.30591061708510
        // nan this is where the squared, unsquared box stuff comes in.
        // makes the density 0.41657 per Dr. Monson's comment in e-mail.
        // defaults.boxSize = 7.018;
        // defaults.boxSize = 100;

        // INTERMOLECULAR POTENTIAL STUFF

        // This potential is the intermolecular potential between atoms on
        // different molecules. We use the class "Potential" because we are
        // reusing the instance as we define each potential.
        Potential potential = new P2HardSphere(space, defaults.atomSize,
                defaults.ignoreOverlap);

        // here, we add the species to the PotentialMaster, using types.
        // The PotentialMaster generates a group potential and automatically
        // does a lot of the stuff which we have to do for the intramolecular
        // potential manually.
        AtomTypeSphere sphereType = (AtomTypeSphere) ((AtomFactoryHomo) species
                .moleculeFactory()).getChildFactory().getType();

        // Add the Potential to the PotentialMaster
        potentialMaster.addPotential(potential,
                new AtomType[] { sphereType, sphereType });

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

        // Initialize the positions of the atoms.
        config.initializeCoordinates(phase);

        integrator.setPhase(phase);

        // nan this will need to be changed
        // pri = new PairIndexerMolecule(phase, new PrimitiveHexane(space));
    }

    public static void main(String[] args) {
        int numMolecules = 144; // 144
        boolean graphic = false;

        // spaces are now singletons; we can only have one instance, so we call
        // it with this method, not a "new" thing.
        TestHexaneCBMCOnly sim = new TestHexaneCBMCOnly(Space3D.getInstance(),
                numMolecules);

        if (graphic) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME);
            simGraphic.makeAndDisplayFrame(APP_NAME);
        } else {
            // PDBWriter write = new PDBWriter(sim.phase);
            // write.setFileName("HexaneCBMCOnly");
            // sim.getController().addAction(write); //after it runs, it writes.

            // IntervalActionAdapter writeAdapter = new
            // IntervalActionAdapter(write, sim.integrator);
            // writeAdapter.setActionInterval(1);
            // sim.integrator.addListener(writeAdapter);
            // writeAdapter.setPriority(0);

            // PrimitiveHexane primitive =
            // (PrimitiveHexane)sim.lattice.getPrimitive();
            // // primitive doesn't need scaling. The boundary was designed to
            // be commensurate with the primitive
            // WaveVectorFactorySimple waveVectorFactory = new
            // WaveVectorFactorySimple(primitive);
            // // we need to set this up now even though we don't use it during
            // equilibration so that
            // // the meter can grab the lattice points
            // MeterNormalMode meterNormalMode = new MeterNormalMode();
            // meterNormalMode.setWaveVectorFactory(waveVectorFactory);
            // meterNormalMode.setCoordinateDefinition(new
            // CoordinateDefinitionHexane());
            // meterNormalMode.setPhase(sim.phase);

            long nSteps = 100;
//            sim.activityIntegrate.setMaxSteps(nSteps / 10);
//            sim.getController().actionPerformed();
//            System.out.println("equilibration finished");

            // ((MCMoveStepTracker)sim.moveMolecule.getTracker()).setTunable(false);
            // ((MCMoveStepTracker)sim.rot.getTracker()).setTunable(false);
            // System.out.println("1");

            sim.getController().reset();
            sim.activityIntegrate.setMaxSteps(nSteps);

//             IntervalActionAdapter checkAdapter = new 
//                  IntervalActionAdapter(check, sim.integrator);
//             checkAdapter.setActionInterval(100);
//             sim.integrator.addListener(checkAdapter);
            //            
            // IntervalActionAdapter writeAdapter = new
            // IntervalActionAdapter(write, sim.integrator);
            // writeAdapter.setActionInterval(100);
            // sim.integrator.addListener(writeAdapter);

            // IntervalActionAdapter adapter = new
            // IntervalActionAdapter(meterNormalMode);
            // adapter.setActionInterval(100);
            // sim.integrator.addListener(adapter);

            sim.getController().actionPerformed();

            // DataGroup normalModeData = (DataGroup)meterNormalMode.getData();
            // normalModeData.TE(1.0/(sim.phase.getSpeciesMaster().moleculeCount()
//                  *meterNormalMode.getCallCount()));
            // int normalDim =
            // meterNormalMode.getCoordinateDefinition().getCoordinateDim();
            //            
            // IVector[] waveVectors = waveVectorFactory.getWaveVectors();
            // double[] coefficients = waveVectorFactory.getCoefficients();
            //            
            // try {
            // FileWriter fileWriterQ = new FileWriter(filename+".Q");
            // FileWriter fileWriterS = new FileWriter(filename+".S");
            // for (int i=0; i<waveVectors.length; i++) {
            // fileWriterQ.write(Double.toString(coefficients[i]));
            // for (int j=0; j<waveVectors[i].getD(); j++) {
            // fileWriterQ.write(" "+waveVectors[i].x(j));
            // }
            // fileWriterQ.write("\n");
            // DataDoubleArray dataS =
            // (DataDoubleArray)normalModeData.getData(i);
            // for (int k=0; k<normalDim; k++) {
            // fileWriterS.write(Double.toString(dataS.getValue(k*normalDim)));
            // for (int l=1; l<normalDim; l++) {
            // fileWriterS.write(" "+dataS.getValue(k*normalDim+l));
            // }
            // fileWriterS.write("\n");
            // }
            // }
            // fileWriterQ.close();
            // fileWriterS.close();
            // }
            // catch (IOException e) {
            // throw new RuntimeException("Oops, failed to write data "+e);
            // }
            // }

        }

    }

    public ActivityIntegrate activityIntegrate;
    public IntegratorMC integrator;
    public Phase phase;
    public BoundaryDeformablePeriodic bdry;
    public CBMCGrowSolidHexane growMolecule;
    public BravaisLattice lattice;


}