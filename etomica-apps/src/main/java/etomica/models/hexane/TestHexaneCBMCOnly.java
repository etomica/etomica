/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * Created on May 24, 2005
 */
package etomica.models.hexane;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.lattice.BravaisLattice;
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
    public ActivityIntegrate activityIntegrate;
    public IntegratorMC integrator;
    public Box box;
    public BoundaryDeformablePeriodic bdry;
    public CBMCGrowSolidHexane growMolecule;
    public BravaisLattice lattice;
    public TestHexaneCBMCOnly(Space _space, int numMolecules) {
        // super(space, false, new PotentialMasterNbr(space, 12.0));
        // super(space, true, new PotentialMasterList(space, 12.0));
        super(_space);
        PotentialMaster potentialMaster = new PotentialMaster();
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
        ConfigurationLattice config = new ConfigurationLattice(lattice, space);


        SpeciesHexane species = new SpeciesHexane(space);
        addSpecies(species);
        bdry = new BoundaryDeformableLattice(primitive, new int[] {4, 6, 6 });
        box = new Box(bdry, space);
        addBox(box);
        box.setNMolecules(species, numMolecules);
        // config.initializeCoordinates(box);

        integrator = new IntegratorMC(potentialMaster, getRandom(), 1.0, box);

        growMolecule = new CBMCGrowSolidHexane(potentialMaster,
                getRandom(), space, integrator, box, species, 20);
        growMolecule.setBox(box);
        integrator.getMoveManager().addMCMove(growMolecule);

        // nan we're going to need some stuff in there to set the step sizes and
        // other stuff like that.

        integrator.setIsothermal(true);
        activityIntegrate = new ActivityIntegrate(integrator);
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
        Potential potential = new P2HardSphere(space);

        // here, we add the species to the PotentialMaster, using types.
        // The PotentialMaster generates a group potential and automatically
        // does a lot of the stuff which we have to do for the intramolecular
        // potential manually.
        AtomType sphereType = species.getLeafType();

        // Add the Potential to the PotentialMaster
        potentialMaster.addPotential(potential,
                new AtomType[]{sphereType, sphereType});

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
        config.initializeCoordinates(box);

        integrator.setBox(box);

        // nan this will need to be changed
        // pri = new PairIndexerMolecule(box, new PrimitiveHexane(space));
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
            // PDBWriter write = new PDBWriter(sim.box);
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
            // meterNormalMode.setBox(sim.box);

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
            // normalModeData.TE(1.0/(sim.box.getSpeciesMaster().moleculeCount()
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


}
