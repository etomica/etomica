package etomica.simulation.prototypes;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.data.meter.MeterEnergy;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntervalActionAdapter;
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
    public DisplayPhase display;
    public DisplayPlot plot;
    public MeterEnergy energy;

    public static void main(String[] args) {
    	LjMd3D sim = new LjMd3D();
    	SimulationGraphic simgraphic = new SimulationGraphic(sim);
    	simgraphic.makeAndDisplayFrame();
    }
    
    public LjMd3D() {
        super(Space3D.getInstance());
        PotentialMaster potentialMaster = new PotentialMaster(this);
        defaults.makeLJDefaults();
        integrator = new IntegratorVelocityVerlet(this, potentialMaster);
        integrator.setTimeStep(0.01);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(this,integrator);
        activityIntegrate.setSleepPeriod(2);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);
        phase.getAgent(species).setNMolecules(50);
        phase = new Phase(this);
        addPhase(phase);
        potential = new P2LennardJones(this);
        potentialMaster.addPotential(potential,new Species[]{species,species});
        
//      elementCoordinator.go();
        //explicit implementation of elementCoordinator activities
        integrator.setPhase(phase);
        PhaseImposePbc imposepbc = new PhaseImposePbc();
        imposepbc.setPhase(phase);
        integrator.addListener(new IntervalActionAdapter(imposepbc));
		
		energy = new MeterEnergy(potentialMaster);
//		energy.setHistorying(true);
//		
//		energy.getHistory().setHistoryLength(500);
//		
//		plot = new DisplayPlot(this);
//		plot.setLabel("Energy");
//		plot.setDataSources(energy.getHistory());
		
    }
    
}