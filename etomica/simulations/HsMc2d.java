package etomica.simulations;
import etomica.*;
import etomica.graphics.*;

/**
 * Simple hard-sphere Monte Carlo simulation in 2D.
 *
 * @author David Kofke
 */
 
public class HsMc2d extends SimulationGraphic {
    
    public IntegratorMC integrator;
    public MCMoveAtom mcMoveAtom;
    public SpeciesSpheresMono species;
    public Phase phase;
    public P2HardSphere potential;
    public Controller controller;
    public DisplayPhase display;
    public MeterCycles meterCycles;
    public DisplayBox displayCycles;

    public HsMc2d() {
 //       super(new Space2D());
        super(new etomica.space.continuum.Space(2));
        Simulation.instance = this;
	    phase = new Phase(this);
	    integrator = new IntegratorMC(this);
	    mcMoveAtom = new MCMoveAtom(integrator);
	    species = new SpeciesSpheresMono(this);
	    potential = new P2HardSphere();
	    controller = new Controller(this);
	    display = new DisplayPhase(this);
	    meterCycles = new MeterCycles(this);
	    displayCycles = new DisplayBox(this,meterCycles);
		panel().setBackground(java.awt.Color.yellow);
		elementCoordinator.go();
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        HsMc2d sim = new HsMc2d();
		sim.elementCoordinator.go(); 
		SimulationGraphic.makeAndDisplayFrame(sim);
	//	sim.controller.start();
    }//end of main
    
}