package etomica.simulations;
import etomica.*;

/**
 * Simple square-well molecular dynamics simulation in 2D
 */
 
public class SwMd2D extends Simulation {
    
    public IntegratorHard integrator;
    public SpeciesDisks species;
    public Phase phase;
    public P2SimpleWrapper p2;
    public PotentialSquareWell potential;
    public Controller controller;
    public DisplayPhase display;

    public SwMd2D() {
        super(new Space2D());
        Simulation.instance = this;
        Default.ATOM_SIZE = 2.0;
	    integrator = new IntegratorHard(this);
	    integrator.setInterval(5);
	    integrator.setSleepPeriod(1);
	    integrator.setTimeStep(0.02);
	    integrator.setTemperature(450.);
	    integrator.setIsothermal(true);
	    species = new SpeciesDisks(this);
	    species.setNMolecules(80);
	    phase = new Phase(this);
	    potential = new PotentialSquareWell(this);
	    p2 = new P2SimpleWrapper(this,potential);
	    controller = new Controller(this);
	    display = new DisplayPhase(this);
	    IntegratorMD.Timer timer = integrator.new Timer(integrator.chronoMeter());
	    timer.setUpdateInterval(10);
		setBackground(java.awt.Color.yellow);
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        javax.swing.JFrame f = new javax.swing.JFrame();   //create a window
        f.setSize(600,350);
        
        Simulation sim = new SwMd2D();
		sim.elementCoordinator.go(); 
		
        f.getContentPane().add(sim.panel());
        
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }//end of main
    
}