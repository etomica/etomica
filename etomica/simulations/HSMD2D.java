package etomica.simulations;
import etomica.*;
import etomica.graphics.*;

/**
 * Simple hard-sphere molecular dynamics simulation in 2D.
 *
 * @author David Kofke
 */
 
public class HSMD2D extends SimulationGraphic {
    
    public IntegratorHard integrator;
    public Species species;
    public Phase phase;
    public Potential2 potential;
    public Controller controller;
    public DisplayPhase display;

    public HSMD2D() {
        super(new etomica.space.continuum.Space(2));
  //      super(new Space2D());
        Simulation.instance = this;
	    integrator = new IntegratorHard(this);
	    species = new SpeciesSpheresMono(this);
	    species.setNMolecules(26);
	    phase = new Phase(this);
	    potential = new P2HardSphere();
	    potential.setSpecies(species,species);
	    controller = new Controller(this);
	    display = new DisplayPhase(this);
	    DisplayTimer timer = new DisplayTimer(integrator);
	    timer.setUpdateInterval(10);
		panel().setBackground(java.awt.Color.yellow);
		elementCoordinator.go();
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        HSMD2D sim = new HSMD2D();
		SimulationGraphic.makeAndDisplayFrame(sim);
	//	sim.controller.start();
    }//end of main
    
}