package etomica.simulations;
import etomica.*;
import etomica.graphics.*;

/**
 * Simple square-well molecular dynamics simulation in 2D
 */
 
public class SwMd2D extends SimulationGraphic {
    
    public IntegratorHard integrator;
    public SpeciesSpheresMono species;
    public Phase phase;
    public P2SquareWell potential;
    public Controller controller;
    public DisplayPhase display;

    public SwMd2D() {
        super(new etomica.space.continuum.Space(2));
        Simulation.instance = this;
        Default.ATOM_SIZE = 2.0;
	    integrator = new IntegratorHard(this);
	    integrator.setInterval(5);
	    integrator.setSleepPeriod(1);
	    integrator.setTimeStep(0.02);
	    integrator.setTemperature(450.);
	    integrator.setIsothermal(true);
	    species = new SpeciesSpheresMono(this);
	    species.setNMolecules(80);
	    phase = new Phase(this);
	    potential = new P2SquareWell();
	    controller = new Controller(this);
	    display = new DisplayPhase(this);
	    DisplayTimer timer = new DisplayTimer(integrator);
	    timer.setUpdateInterval(10);
		panel().setBackground(java.awt.Color.yellow);
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        SwMd2D sim = new SwMd2D();
		sim.elementCoordinator.go();
		sim.makeAndDisplayFrame();
    }//end of main
    
}