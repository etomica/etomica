package etomica.simulation.prototypes;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.DisplayBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.phase.Phase;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
 
public class LjMd3D extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public IntegratorVelocityVerlet integrator;
    public SpeciesSpheresMono species;
    public Phase phase;
    public P2LennardJones potential;
    public Controller controller;
    public MeterPotentialEnergy energy;
    public AccumulatorAverage avgEnergy;

    public static void main(String[] args) {
    	final String APP_NAME = "LjMd3D";
    	final LjMd3D sim = new LjMd3D();
    	final SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME, 3);

        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getDisplayPhasePaintAction(sim.phase));

        simGraphic.makeAndDisplayFrame(APP_NAME);

        DisplayBoxesCAE display = new DisplayBoxesCAE();
        display.setAccumulator(sim.avgEnergy);
        simGraphic.add(display);
    }
    
    public LjMd3D() {
        super(Space3D.getInstance());
        PotentialMaster potentialMaster = new PotentialMaster(space);
        double sigma = 3.0;
        integrator = new IntegratorVelocityVerlet(this, potentialMaster);
        integrator.setTimeStep(0.01);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator, true, false);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);
        phase = new Phase(this);
        addPhase(phase);
        phase.getAgent(species).setNMolecules(50);
        potential = new P2LennardJones(space, sigma, 1.0);
        potentialMaster.addPotential(potential,new Species[]{species,species});
        
        integrator.setPhase(phase);
        PhaseImposePbc imposepbc = new PhaseImposePbc();
        imposepbc.setPhase(phase);
        integrator.addIntervalAction(imposepbc);
		
        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc());
        configuration.initializeCoordinates(phase);
        energy = new MeterPotentialEnergy(potentialMaster);
        energy.setPhase(phase);
        avgEnergy = new AccumulatorAverage(200);
        avgEnergy.setPushInterval(10);
        DataPump pump = new DataPump(energy, avgEnergy);
        integrator.addIntervalAction(pump);
        integrator.setActionInterval(pump, 10);
        register(energy, pump);
    }
    
}