package etomica.simulations;
import etomica.*;

/**
 * Simple hard-sphere Monte Carlo simulation in 2D
 */
 
public class HsMc2d extends Simulation {
    
    public IntegratorMC integrator;
    public MCMoveAtom mcMoveAtom;
    public SpeciesDisks species;
    public Phase phase;
    public P2SimpleWrapper p2;
    public PotentialHardDisk potential;
    public Controller controller;
    public DisplayPhase display;
    public MeterCycles meterCycles;
    public DisplayBox displayCycles;

    public HsMc2d() {
        super(new Space2D());
        Simulation.instance = this;
	    integrator = new IntegratorMC(this);
	    mcMoveAtom = new MCMoveAtom();
	    integrator.add(mcMoveAtom);
	    species = new SpeciesDisks(this);
	    phase = new Phase(this);
	    potential = new PotentialHardDisk(this);
	    p2 = new P2SimpleWrapper(this,potential);
	    controller = new Controller(this);
	    display = new DisplayPhase(this);
	    meterCycles = new MeterCycles(this);
	    displayCycles = new DisplayBox(this,meterCycles);
		setBackground(java.awt.Color.yellow);
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        javax.swing.JFrame f = new javax.swing.JFrame();   //create a window
        f.setSize(600,350);
        
        Simulation sim = new HsMc2d();
		sim.elementCoordinator.go(); 
		
        f.getContentPane().add(sim.panel());
        
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }//end of main
    
}