package etomica.simulations;
import etomica.Controller;
import etomica.Default;
import etomica.Phase;
import etomica.Simulation;
import etomica.Species;
import etomica.SpeciesSpheresMono;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.data.meter.MeterEnergy;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.P2LennardJones;
import etomica.space3d.Space3D;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
 
public class LjMd3D extends Simulation {
    
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
        Default.makeLJDefaults();
        integrator = new IntegratorVelocityVerlet(potentialMaster, space);
        integrator.setTimeStep(0.01);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(2);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
        species.setNMolecules(50);
        phase = new Phase(this);
        potential = new P2LennardJones(space);
        this.potentialMaster.setSpecies(potential,new Species[]{species,species});
        
//      elementCoordinator.go();
        //explicit implementation of elementCoordinator activities
        integrator.addPhase(phase);
        PhaseImposePbc imposepbc = new PhaseImposePbc();
        imposepbc.setPhase(phase);
        integrator.addListener(imposepbc);
		
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