package etomica.simulation.prototypes;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.config.ConfigurationSequential;
import etomica.graphics.DisplayPhase;
import etomica.integrator.IntegratorHard;
import etomica.phase.Phase;
import etomica.potential.P2SquareWell;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.util.Default;

/**
 * Simple square-well molecular dynamics simulation in 2D
 */
 
public class SwMd2D extends Simulation {
    
    public IntegratorHard integrator;
    public SpeciesSpheresMono species;
    public Phase phase;
    public P2SquareWell potential;
    public Controller controller;
    public DisplayPhase display;

    public SwMd2D() {
        this(Space2D.getInstance());
    }
    public SwMd2D(Space space) {
        super(space);
        Default.makeLJDefaults();
        Default.ATOM_SIZE = 0.8;
        integrator = new IntegratorHard(potentialMaster);
        integrator.setTimeStep(0.01);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setInterval(5);
        activityIntegrate.setSleepPeriod(1);
        integrator.setTimeStep(0.02);
        integrator.setTemperature(300.);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
        species.setNMolecules(50);
        phase = new Phase(this);
        new ConfigurationSequential(space).initializeCoordinates(phase);
        potential = new P2SquareWell(space);
        this.potentialMaster.setSpecies(potential,new Species[]{species,species});
        
//      elementCoordinator.go();
        //explicit implementation of elementCoordinator activities
        integrator.addPhase(phase);
        integrator.addListener(new PhaseImposePbc(phase));
    } 
}