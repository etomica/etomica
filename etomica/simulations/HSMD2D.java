package etomica.simulations;
import etomica.*;

/**
 * Simple hard-sphere molecular dynamics simulation in 2D
 */
 
public class HSMD2D extends Simulation {
    
    public IntegratorHard integrator;
    public Species species;
    public Phase phase;
    public Potential2 potential;
    public Controller controller;
    public DisplayPhase display;

    public HSMD2D() {
        super(new Space2D());
        Simulation.instance = this;
	    integrator = new IntegratorHard(this);
	    species = new SpeciesSpheresMono(this);
	    species.setNMolecules(25);
	    phase = new Phase(this);
	    potential = new P2HardSphere(this);
	    controller = new Controller(this);
	    display = new DisplayPhase(this);
	    IntegratorMD.Timer timer = integrator.new Timer(integrator.chronoMeter());
	    timer.setUpdateInterval(10);
		setBackground(java.awt.Color.yellow);
		elementCoordinator.go();
		
        potential.setIterator(new AtomPairIterator(phase));
        potential.set(species.getAgent(phase));
		
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        HSMD2D sim = new HSMD2D();
		Simulation.makeAndDisplayFrame(sim);
	//	sim.controller.start();
    }//end of main
    
}