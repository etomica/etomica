package etomica.simulations;
import etomica.*;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 2D
 */
 
public class LjMd2D extends Simulation {
    
    public IntegratorVelocityVerlet integrator;
    public SpeciesDisks species;
    public Phase phase;
    public P2SimpleWrapper p2;
    public PotentialLJ potential;
    public Controller controller;
    public DisplayPhase display;
    public DisplayPlot plot;

    public LjMd2D() {
        super(new Space2D());
        Simulation.instance = this;
	    integrator = new IntegratorVelocityVerlet(this);
	    species = new SpeciesDisks(this);
	    phase = new Phase(this);
	    potential = new PotentialLJ();
	    p2 = new P2SimpleWrapper(this,potential);
	    controller = new Controller(this);
	    display = new DisplayPhase(this);
	    IntegratorMD.Timer timer = integrator.new Timer(integrator.chronoMeter());
	    timer.setUpdateInterval(10);
		setBackground(java.awt.Color.yellow);
		
		phase.energy.setHistorying(true);
		phase.energy.setActive(true);
		
		phase.energy.getHistory().setNValues(500);
		
		plot = new DisplayPlot(this);
		plot.setLabel("Energy");
		plot.setDataSource(phase.energy.getHistory());
		
		integrator.setSleepPeriod(2);
		
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        javax.swing.JFrame f = new javax.swing.JFrame();   //create a window
        f.setSize(600,350);
        
        Simulation sim = new LjMd2D();
		sim.elementCoordinator.go(); 
		
        f.getContentPane().add(sim.panel());
        
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }//end of main
    
}