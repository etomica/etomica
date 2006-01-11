/*
 * Created on May 24, 2005
 */
package etomica.models.hexane;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.atom.iterator.ApiBuilder;
import etomica.atom.iterator.ApiIntragroup;
import etomica.atom.iterator.AtomIteratorArrayList;
import etomica.atom.iterator.IteratorDirective;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressureHard;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayPhase;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntervalActionAdapter;
import etomica.nbr.CriterionBondedSimple;
import etomica.nbr.CriterionMolecular;
import etomica.nbr.CriterionMolecularNonAdjacent;
import etomica.nbr.CriterionSimple;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.phase.Phase;
import etomica.potential.P2HardBond;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential;
import etomica.potential.Potential2;
import etomica.potential.PotentialGroup;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.util.Default;

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
    public IntegratorHard integrator;

    public Phase phase;

    public BoundaryDeformablePeriodic bdry;

    public TestHexane(Space space, int numMolecules) {
        //super(space, false, new PotentialMasterNbr(space, 12.0));
        super(space, true, new PotentialMasterList(space, 12.0));
        int chainLength = 6;
        int numAtoms = numMolecules * chainLength;
        ConfigurationHexane config = new ConfigurationHexane(space);

        //This is the factor that multiples by the range of the potential in
        // order to define the area/volume in which neighbors are searched for.
        //This becomes the bond delta, which is the percentage the bond can
        // stretch, and I assume compress.
        double neighborRangeFac = 1.2;

        //The bondfactor is not used, so it is set to 0.
        double bondFactor = 0.0;
        Default def = new Default();
        def.makeLJDefaults();

        double timeStep = 0.005;
        double simTime = 100000.0 / numAtoms;
        int nSteps = (int) (simTime / timeStep);

        //nan The box size we want is 5.72906360610622 by 11.21417818673970 by
        // 7.30591061708510
        //nan this is where the squared, unsquared box stuff comes in.
        //makes the density 0.41657 per Dr. Monson's comment in e-mail.
        def.boxSize = 7.018;
        integrator = new IntegratorHard(potentialMaster, timeStep,
                def.temperature);
        integrator.setTimeStep(timeStep);
        integrator.setIsothermal(true);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(this,
                integrator);
        //Added to make it work.
        NeighborListManager nbrManager = ((PotentialMasterList) potentialMaster)
                .getNeighborManager();
        integrator.addListener(nbrManager);
        nbrManager.setRange(defaults.atomSize * neighborRangeFac);
        ((PotentialMasterList) potentialMaster).setCellRange(2);
        ((PotentialMasterList) potentialMaster).setRange(neighborRangeFac
                * defaults.atomSize);

        getController().addAction(activityIntegrate);
        activityIntegrate.setMaxSteps(nSteps);

        //INTERMOLECULAR POTENTIAL STUFF

        //This potential is the intermolecular potential between atoms on
        // different molecules. We use the class "Potential" because we are
        // reusing the instance as we define each potential.
        Potential potential = new P2HardSphere(space, def.atomSize,
                def.ignoreOverlap);

        //here, we add the species to the PotentialMaster, using types.
        //The PotentialMaster generates a group potential and automatically
        // does a lot of the stuff which we have to do for the intramolecular
        // potential manually.
        SpeciesHexane species = new SpeciesHexane(this);
        species.setNMolecules(numMolecules);
        AtomTypeSphere sphereType = (AtomTypeSphere) ((AtomFactoryHomo) species
                .moleculeFactory()).getChildFactory().getType();

        //Using a NeighborCriterion allows us to re-use the handle nbrCriterion
        // for other types
        //The NeighborCriterion defines the neighbors on the neighbor list.
        //It is a range in the potential, multiplied by some factor to extend
        // it (neighborRangeFac).
        NeighborCriterion nbrCriterion = new CriterionSimple(this, potential
                .getRange(), neighborRangeFac * potential.getRange());

        //here we create an on-the-molecule-checking sub-criterion/criterion
        // wrapper.
        CriterionMolecular sameMoleculeCrit = new CriterionMolecular(
                nbrCriterion);
        sameMoleculeCrit.setIntraMolecular(false);

        //we add the new criterion to the potential.
        ((Potential2) potential).setCriterion(sameMoleculeCrit);
        nbrManager.addCriterion(sameMoleculeCrit, new AtomType[] {sphereType,sphereType});

        //Add the Potential to the PotentialMaster
        potentialMaster.addPotential(potential, new AtomType[] { sphereType,
                sphereType });
        
         //INTRAMOLECULAR POTENTIAL STUFF

        //This PotentialGroup will hold all the intramolecular potentials.
        //We give 1 as the argument because we are using 1 molecule to iterate
        // on. The actual interactions between the atoms on the molecules will
        // be calculated by a Potential2, but their summation is the molecule's
        //effect on itself, which is a Potential1, or a Potential with nBody =
        // 1.
        PotentialGroup potentialChainIntra = new PotentialGroup(1, space);

        //BONDED INTERACTIONS

        //nan This potential simulates the bonds between atoms in a molecule.
        // It will be superceded by a set of MC moves at some point in the
        // future.
        //This potential uses hard sphere interactions to model the bonded
        // interactions
        // of the atoms of the molecule.
        potential = new P2HardBond(space, def.atomSize, neighborRangeFac,
                def.ignoreOverlap);

        //Only the atoms next to each other interact, so we have two criteria:
        //		The atoms must be on the same molecule- CriterionMolecular
        //		The atoms must be neighbors
        //Both are included in the CriterionBondedSimple class, so we use it.
        NeighborCriterion criterion = new CriterionBondedSimple(nbrCriterion);
        ((Potential2) potential).setCriterion(criterion);
        

        //We will need an atom pair iterator (Api) that runs through the atoms
        // on a single molecule.
        //The atom pair iterator (Api) runs through the atoms on a single
        // molecule.
        //  It has an inner loop and an outer loop.
        ApiIntragroup bonded = ApiBuilder.makeAdjacentPairIterator();

        //We add the Potential and its Iterator to the PotentialGroup, in one
        // fell swoop. Yay us!
        potentialChainIntra.addPotential(potential, bonded);

        //NONBONDED INTERACTIONS
        //This potential describes the basic hard sphere interactions between
        // 2 atoms of a molecule.
        potential = new P2HardSphere(space, def.atomSize, false);

        //Only the atoms next to each other interact, so we have two criteria:
        //		The atoms must be on the same molecule- CriterionMolecular
        //		The atoms must be separated by 3 bonds.
        //We end up needing to do stuff with the
        criterion = new CriterionMolecularNonAdjacent(3, nbrCriterion);
        ((Potential2) potential).setCriterion(criterion);

        //This iterator runs through the atoms on a single molecule. It skips
        // the specified number of atoms.
        //The inner iterator runs through
        //The outer iterator runs through
        //Similar to ApiBuilder.makeNonAdjacentPairIterator
        AtomIteratorArrayList aiInnerUp = new AtomIteratorArrayList(IteratorDirective.UP, 3);
        AtomIteratorArrayList aiInnerDn = new AtomIteratorArrayList(IteratorDirective.DOWN, 3);
        ApiIntragroup intra = new ApiIntragroup(aiInnerUp, aiInnerDn);

        //Now we add this potential to the PotentialGroup.
        potentialChainIntra.addPotential(potential, intra);

        //Now we add this PotentialGroup to the PotentialMaster.
        potentialMaster.addPotential(potentialChainIntra,
                new AtomType[] { species.getMoleculeType() });

        //nan, um, how do we worry about running through all the molecules?

        phase = new Phase(this);
        //nan we need to do more with the boundary
        bdry =  new BoundaryHexane(space);
        phase.setBoundary(bdry);

        integrator.setPhase(phase);

        //Initialize the positions of the atoms.
        config.initializeCoordinates(phase);

    }

    public static void main(String[] args) {
        int numMolecules = 144;

        //spaces are now singletons; we can only have one instance, so we call
        // it with this method, not a "new" thing.
        TestHexane sim = new TestHexane(Space3D.getInstance(), numMolecules);

        System.out.println("Happy Goodness!!");

        MeterPressureHard pMeter = new MeterPressureHard(sim.space,
                sim.integrator);
        pMeter.setPhase(sim.phase);
        MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(
                sim.integrator);
        energyMeter.setPhase(sim.phase);
        AccumulatorAverage energyAccumulator = new AccumulatorAverage(sim);
        DataPump energyManager = new DataPump(energyMeter, energyAccumulator);
        energyAccumulator.setBlockSize(50);
        pMeter.setPhase(sim.phase);
        //This is the wire to plug the pump into the electric socket and make
        // it run.
        new IntervalActionAdapter(energyManager, sim.integrator);

        //This starts the simulation, sets up its own thread to start the sim,
        // and then moves on down the list.
        //sim.getController().actionPerformed();

        double avgPE = ((DataDouble) ((DataGroup) energyAccumulator.getData())
                .getData(AccumulatorAverage.AVERAGE.index)).x;
        System.out.println("PE  " + avgPE);
        
        SimulationGraphic simGraphic = new SimulationGraphic(sim);
        simGraphic.makeAndDisplayFrame();

    }

}