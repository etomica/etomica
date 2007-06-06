package etomica.models.hexane;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntervalActionAdapter;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.BravaisLattice;
import etomica.lattice.crystal.Primitive;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.MCMoveMoleculeCoupled;
import etomica.normalmode.MeterNormalMode;
import etomica.normalmode.WaveVectorFactory;
import etomica.normalmode.WaveVectorFactorySimple;
import etomica.normalmode.WriteS;
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

public class TestHexaneLowMedium extends Simulation {

    public ActivityIntegrate activityIntegrate;
    public IntegratorMC integrator;
    public Phase phase;
    public BoundaryDeformablePeriodic bdry;
    public BravaisLattice lattice;
    public CoordinateDefinition coordinateDefinition;
    public Primitive primitive;

    public MCMoveRotateMolecule3D rot;
    public MCMoveMoleculeCoupled coupledMove;
    public MCMoveCombinedCbmcTranslation cctMove;
    protected CBMCGrowSolidHexane growMolecule;
    private double density;
    
    public TestHexaneLowMedium(Space space, int numMolecules) {
        super(space, false);
        PotentialMaster potentialMaster = new PotentialMaster(this);
        int chainLength = 6;
        int numAtoms = numMolecules * chainLength;
        primitive = new PrimitiveHexane(space);
        // close packed density is 0.4165783882178116
        // Monson reports data for 0.373773507616 and 0.389566754417
        density = 0.373773507616;
        primitive.scaleSize(Math.pow(0.4165783882178116/density,1.0/3.0));
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
        int[] nCells = new int[]{4,3,3};
        bdry = new BoundaryDeformableLattice(primitive, getRandom(), nCells);
        phase = new Phase(bdry);
        addPhase(phase);
        phase.getAgent(species).setNMolecules(numMolecules);
//            config.initializeCoordinates(phase);
        integrator = new IntegratorMC(potentialMaster, getRandom(),
                defaults.temperature);
        
        rot = new MCMoveRotateMolecule3D(potentialMaster, getRandom());
        rot.setPhase(phase);
        rot.setStepSize(0.042);
        integrator.getMoveManager().addMCMove(rot);
        ((MCMoveStepTracker)rot.getTracker()).setNoisyAdjustment(true);
        
        coupledMove = new MCMoveMoleculeCoupled(potentialMaster, getRandom());
        integrator.getMoveManager().addMCMove(coupledMove);
        
        growMolecule = new CBMCGrowSolidHexane(potentialMaster,
                getRandom(), integrator, phase, species, 20);
        growMolecule.setPhase(phase);
        cctMove = new MCMoveCombinedCbmcTranslation(potentialMaster, growMolecule, getRandom());
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
//            defaults.boxSize = 7.018;
//            defaults.boxSize = 100;

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


        //Initialize the positions of the atoms.
        coordinateDefinition = new CoordinateDefinitionHexane(phase, primitive, species);
        coordinateDefinition.initializeCoordinates(nCells);

        integrator.setPhase(phase);
        
    }

    public static void main(String[] args) {
        //defaults
        int D = 3;
        int nA = 36;
        double density = 0.37;
        long nSteps = 100;
        int numMolecules = nA; //144
        
        String filename = "normal_modes"+nA+"_"+((int)(density*100))+"hexane"+nSteps;

        System.out.println("Running hard sphere hexane simulation");
        System.out.println(numMolecules + " molecules at density " + density);
        System.out.println("output data to " + filename);

        //spaces are now singletons; we can only have one instance, so we call
        // it with this method, not a "new" thing.
        TestHexaneLowMedium sim = new TestHexaneLowMedium(Space3D.getInstance(), numMolecules);

        PrimitiveHexane primitive = (PrimitiveHexane)sim.lattice.getPrimitive();
        
        // set up normal-mode meter
        MeterNormalMode meterNormalMode = new MeterNormalMode();
        meterNormalMode.setCoordinateDefinition(sim.coordinateDefinition);
        WaveVectorFactory waveVectorFactory;
        waveVectorFactory = new WaveVectorFactorySimple(primitive);
        
        meterNormalMode.setWaveVectorFactory(waveVectorFactory);
        meterNormalMode.setPhase(sim.phase);

        IntervalActionAdapter fooAdapter = new IntervalActionAdapter(
                meterNormalMode);
        fooAdapter.setActionInterval(numMolecules);
        sim.integrator.addListener(fooAdapter);
        
        sim.activityIntegrate.setMaxSteps(nSteps/10);
        sim.getController().actionPerformed();
        System.out.println("equilibration finished");
        meterNormalMode.reset();
        sim.getController().reset();
        
        sim.activityIntegrate.setMaxSteps(nSteps);      
        sim.getController().actionPerformed();
 
        System.out.println("Calculating all the fun eigenstuff");
        
        WriteS sWriter = new WriteS();
        sWriter.setFilename(filename);
        sWriter.setMeter(meterNormalMode);
        sWriter.setWaveVectorFactory(waveVectorFactory);
        sWriter.setOverwrite(true);
        sWriter.actionPerformed();
        
        System.out.println(filename + " is finished running");
        System.out.println(nSteps + "  steps");
    }

}


