package etomica.simulations;
import etomica.*;
import etomica.graphics.*;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 2D
 */
 
public class Lj2D_NoIntegrator extends SimulationGraphic {
    
    public SpeciesSpheresMono species;
    public Phase phase;
    public P2LennardJones potential;
    public Controller controller;
    public DisplayPhase display;
    public DisplayPlot plot;
    public MeterEnergy energy;

    public Lj2D_NoIntegrator() {
        super(new Space2D());
        Simulation.instance = this;
	    species = new SpeciesSpheresMono(this);
	    phase = new Phase(this);
	    potential = new P2LennardJones();
	    controller = new Controller(this);
	    display = new DisplayPhase(this);
		panel().setBackground(java.awt.Color.yellow);
				
		energy = new MeterEnergy(this);
		energy.setHistorying(true);
		energy.setActive(true);
		
		energy.getHistory().setHistoryLength(500);

		plot = new DisplayPlot(this);
		plot.setLabel("Energy");
		plot.setDataSource(energy.getHistory());
				
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        SimulationGraphic sim = new Lj2D_NoIntegrator();
        new IntegratorGear4(sim);
		sim.elementCoordinator.go(); 
		
        sim.makeAndDisplayFrame();
    }//end of main
    
}