package etomica.simulations;
import etomica.Controller;
import etomica.Default;
import etomica.Phase;
import etomica.Simulation;
import etomica.Space2D;
import etomica.Species;
import etomica.SpeciesSpheresMono;
import etomica.action.activity.ActivityIntegrate;
import etomica.data.meter.MeterEnergy;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayPlot;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.P2LennardJones;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 2D
 */
 
public class LjMd2D extends Simulation {
    
    public IntegratorVelocityVerlet integrator;
    public SpeciesSpheresMono species;
    public Phase phase;
    public P2LennardJones potential;
    public Controller controller;
    public DisplayPhase display;
    public DisplayPlot plot;
    public MeterEnergy energy;

    public LjMd2D() {
        super(Space2D.INSTANCE);
        Default.makeLJDefaults();
        integrator = new IntegratorVelocityVerlet(potentialMaster, space);
        integrator.setTimeStep(0.01);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(2);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
        species.setNMolecules(50);
        phase = new Phase(space);
        potential = new P2LennardJones(space);
        this.potentialMaster.setSpecies(potential,new Species[]{species,species});
        
//      elementCoordinator.go();
        //explicit implementation of elementCoordinator activities
        phase.speciesMaster.addSpecies(species);
        integrator.addPhase(phase);
		
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