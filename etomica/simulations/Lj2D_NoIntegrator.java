package etomica.simulations;
import etomica.*;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 2D
 */
 
public class Lj2D_NoIntegrator extends Simulation {
    
    public SpeciesDisks species;
    public Phase phase;
    public P2SimpleWrapper p2;
    public PotentialLJ potential;
    public Controller controller;
    public DisplayPhase display;
    public DisplayPlot plot;

    public Lj2D_NoIntegrator() {
        super(new Space2D());
        Simulation.instance = this;
	    species = new SpeciesDisks(this);
	    phase = new Phase(this);
	    potential = new PotentialLJ();
	    p2 = new P2SimpleWrapper(this,potential);
	    controller = new Controller(this);
	    display = new DisplayPhase(this);
		setBackground(java.awt.Color.yellow);
		
		phase.energy.setHistorying(true);
		phase.energy.setActive(true);
		
		phase.energy.getHistory().setNValues(500);
		
		plot = new DisplayPlot(this);
		plot.setLabel("Energy");
		plot.setDataSource(phase.energy.getHistory());
				
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        javax.swing.JFrame f = new javax.swing.JFrame();   //create a window
        f.setSize(600,350);
        
        Simulation sim = new Lj2D_NoIntegrator();
        new IntegratorGear4(sim);
		sim.elementCoordinator.go(); 
		
        f.getContentPane().add(sim.panel());
        
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }//end of main
    
}