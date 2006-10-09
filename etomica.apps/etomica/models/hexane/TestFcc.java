package etomica.models.hexane;

import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.DataLogger;
import etomica.data.DataPump;
import etomica.data.DataTableWriter;
import etomica.data.types.CastGroupToDoubleArray;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntervalActionAdapter;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.Primitive;
import etomica.phase.Phase;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;

/**
 * A monatomic fcc hard sphere simulation to test a new energy method.
 * 
 * @author nancycribbin
 * 
 */
public class TestFcc extends Simulation {

    public TestFcc(Space space, int numAtoms) {
        super(space, true, new PotentialMaster(space));

        LatticeCubicFcc lattice = new LatticeCubicFcc(1.0);
        ConfigurationLattice config = new ConfigurationLattice(lattice);
//        config.setRescalingToFitVolume(false);
        config.setRememberingIndices(true);

        defaults.makeLJDefaults();
        defaults.atomSize = 1.0;
        defaults.ignoreOverlap = true;

        SpeciesSpheresMono species = new SpeciesSpheresMono(this);

        phase = new Phase(this);
        phase.getAgent(species).setNMolecules(numAtoms);

        config.initializeCoordinates(phase);

        integrator = new IntegratorHard(this);

        integrator.setIsothermal(false);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(this,
                integrator);
        activityIntegrate.setMaxSteps(2000000);
        getController().addAction(activityIntegrate);

        double timeStep = 0.005;
        double simTime = 100.0 / numAtoms;
        int nSteps = (int) (simTime / timeStep);

        getController().addAction(activityIntegrate);
        // activityIntegrate.setMaxSteps(nSteps);

        Potential potential = new P2HardSphere(space, defaults.atomSize, false);
        AtomTypeSphere sphereType = (AtomTypeSphere) ((AtomFactoryMono) species
                .moleculeFactory()).getType();
        potentialMaster.addPotential(potential, new AtomType[] { sphereType,
                sphereType });

        bdry = new BoundaryRectangularPeriodic(this);
        phase.setBoundary(bdry);

        PhaseImposePbc makeperiodic = new PhaseImposePbc(phase);
        integrator.addListener(makeperiodic);

        config.initializeCoordinates(phase);
        
        //nan this section is a patch
        //first we find out the scaling used in ConfigurationLattice/LatticeCubicFcc
        //then, we create a primitive fcc lattice, and scale it so we can use it in pri.
        ConfigurationLattice.MyLattice myLattice = (ConfigurationLattice.MyLattice)config.getLatticeMemento();
        Vector scaling = myLattice.latticeScaling;
        scaling.TE(0.5*Math.sqrt(2.0)); //we need this because the fccPrimitive.setSize method uses the unit vectors, not the actual lattice vectors, and this eliminates the radical 2 in the unit vectors.
        Primitive fccPrimitive = lattice.getPrimitiveFcc();
        fccPrimitive.setSize(scaling.toArray());
      //nan  phase.setDensity(1.04);
        integrator.setPhase(phase);

        pri = new PairIndexerMolecule(phase, fccPrimitive);
        DataProcessorArrayFlatten squisher = new DataProcessorArrayFlatten();
        meterCorrelation = new MeterCorrelationMatrix(phase, pri);
        AccumulatorAverage accumulator = new AccumulatorAverage(this);
        
        DataPump pump = new DataPump(meterCorrelation, accumulator);
        IntervalActionAdapter adapter = new IntervalActionAdapter(pump, integrator);
        CastGroupToDoubleArray mormon = new CastGroupToDoubleArray();
        DataLogger logger = new DataLogger();
        
        logger.setDataSink(new DataTableWriter());
        logger.setFileName("Happy.txt");
        logger.setAppending(false);
        
        accumulator.addDataSink(mormon, 
                new AccumulatorAverage.StatType[] {AccumulatorAverage.StatType.AVERAGE});
      
        mormon.setDataSink(squisher);
        squisher.setDataSink(logger);
        
        
        System.out.println(phase.getDensity());
        phase.setDensity(1.04);
    }

    /**
     * @param args
     */
    public static void main(String[] args) {
        int nA = 108;
        TestFcc sim = new TestFcc(Space3D.getInstance(), nA);

        SimulationGraphic simG = new SimulationGraphic(sim);
        simG.makeAndDisplayFrame();

//        MeterNormalMode mnm = new MeterNormalMode();
        
        
        
    }

    public IntegratorMD integrator;
    public MeterCorrelationMatrix meterCorrelation;
    public Phase phase;

    public BoundaryRectangularPeriodic bdry;

    public PairIndexerMolecule pri;

}
