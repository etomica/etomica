package etomica.simulation.prototypes;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeSphere;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataPump;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntervalActionAdapter;
import etomica.lattice.Crystal;
import etomica.lattice.LatticeCrystal;
import etomica.lattice.crystal.BasisBetaSnA5;
import etomica.lattice.crystal.PrimitiveTetragonal;
import etomica.phase.Phase;
import etomica.potential.EmbeddedAtomMethodP2;
import etomica.potential.EmbeddedAtomMethodPInitial;
import etomica.potential.EmbeddedAtomMethodPMany;
import etomica.potential.ParameterSetEAM;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.units.ElectronVolt;

/**
 * Molecular-Dynamics Simulation Using the Embedded-Atom Method (EAM) Potential.  
 * 
 * The EAM potential is intended for use with metallic systems.
 * 
 * Because the EAM potential consists of a pair-wise term and a many-body term, it 
 * was necessary to calculate these parts of the potential separately.  The many-body
 * term (the embedding energy function) consists of pair-wise sums multiplied together.
 * These pair-wise terms are calculated in the same class as the pair-wise potential.
 * Both are calculated using the infrastructure already in place for looping through pairs.
 * An allatomAgents array is used to hold the pair-wise terms for the many-body 
 * potential so that the many-body potential may be calculated by referring only to 
 * each atom's information in the array (so that the many-body potential part may be
 * calculated as a one-body potential).
 * 
 * This class was adapted from LjMd3D.java by K.R. Schadel and A. Schultz July 2005.
 */
 
public class EAMMd3D extends Simulation {
    
    public IntegratorVelocityVerlet integrator;
    public SpeciesSpheresMono species;
    public Phase phase;
    public EmbeddedAtomMethodPInitial potential0;
    public EmbeddedAtomMethodP2 potentialA;
    public EmbeddedAtomMethodPMany potentialB;
    public Controller controller;
    public DisplayPhase display;
    public DisplayPlot plot;
    public MeterEnergy energy;
    public ActivityIntegrate activityIntegrate;
    public DataInfo info2;

    public static void main(String[] args) {
    	EAMMd3D sim = new EAMMd3D();
    	MeterPotentialEnergy energyMeter = new MeterPotentialEnergy(sim.potentialMaster);
    	energyMeter.setPhase(sim.phase);
    	AccumulatorAverage energyAccumulator = new AccumulatorAverage(sim);
        DataPump energyManager = new DataPump(energyMeter,energyAccumulator);
        energyAccumulator.setBlockSize(50);
        IntervalActionAdapter adapter = new IntervalActionAdapter(energyManager, sim.integrator);
        adapter.setActionInterval(5);
        SimulationGraphic simgraphic = new SimulationGraphic(sim);
    	simgraphic.makeAndDisplayFrame();
    	ColorSchemeByType.setColor(sim.species.getFactory().getType(), java.awt.Color.red);
    	sim.activityIntegrate.setMaxSteps(100);
    	sim.getController().run();
    	Data data = energyAccumulator.getData(); // kmb change type to Data instead of double[]
        double PE = AccumulatorAverage.AVERAGE.index/sim.species.getAgent(sim.phase).getNMolecules();  // kmb changed 8/3/05
        //double PE = data[AccumulatorAverage.AVERAGE.index]/sim.species.getAgent(sim.phase).getNMolecules();  // orig line
        System.out.println("PE/epsilon="+ElectronVolt.UNIT.fromSim(PE));
    }
    
    public EAMMd3D() {
        super(Space3D.getInstance()); //INSTANCE); kmb change 8/3/05
        integrator = new IntegratorVelocityVerlet(this);
        integrator.setTimeStep(0.01);
        activityIntegrate = new ActivityIntegrate(this,integrator);
        activityIntegrate.setSleepPeriod(2);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
        species.setNMolecules(216);
        ((AtomTypeLeaf)species.getFactory().getType()).setMass(118.71);
        //The diameter given is really the equilibrium distance between the atoms.
        ((AtomTypeSphere)species.getFactory().getType()).setDiameter(3.44); 
        // This forces the EmbeddedAtomMethodP2 to 
        //request an agentIndex from Atom.
        potentialA = new EmbeddedAtomMethodP2(space, ParameterSetEAM.Sn);
        phase = new Phase(this);
        PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.8318, 3.1819);
        LatticeCrystal crystal = new LatticeCrystal(new Crystal(
        		primitive, new BasisBetaSnA5(primitive)));
        Configuration config = new ConfigurationLattice(crystal);
//        phase.setConfiguration(config);  // kmb remove 8/3/05
        config.initializeCoordinates(phase);  // kmb added 8/3/05
        potential0 = new EmbeddedAtomMethodPInitial(space, potentialA);
        potentialB = new EmbeddedAtomMethodPMany(space, ParameterSetEAM.Sn, potentialA);
        this.potentialMaster.setSpecies(potentialB, new Species[]{species});
        this.potentialMaster.setSpecies(potentialA, new Species[]{species,species});
        this.potentialMaster.setSpecies(potential0, new Species[]{species});    

        integrator.setPhase(phase);
        PhaseImposePbc imposepbc = new PhaseImposePbc();
        imposepbc.setPhase(phase);
        integrator.addListener(imposepbc);
		
		energy = new MeterEnergy(potentialMaster);
    }
    
}