package etomica.simulations;
import etomica.Controller;
import etomica.Default;
import etomica.Phase;
import etomica.Simulation;
import etomica.Space;
import etomica.Species;
import etomica.SpeciesSpheresMono;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.graphics.DisplayPhase;
import etomica.integrator.IntegratorHard;
import etomica.potential.P2SquareWell;
import etomica.space2d.Space2D;

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
        this(Space2D.INSTANCE);
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
        potential = new P2SquareWell(space);
        this.potentialMaster.setSpecies(potential,new Species[]{species,species});
        
//      elementCoordinator.go();
        //explicit implementation of elementCoordinator activities
        integrator.addPhase(phase);
        integrator.addIntervalListener(new PhaseImposePbc(phase));
    } 
}