package simulate.simulations;
import simulate.*;

/**
 * Simple hard-sphere molecular dynamics simulation in 2D
 */
 
public class HSMD2D extends Simulation {
    
    public HSMD2D() {
        super(new Space2D());
        Simulation.instance = this;
	    IntegratorHard integratorHard1 = new IntegratorHard(this);
	    SpeciesDisks speciesDisks1 = new SpeciesDisks(this);
	    Phase phase1 = new Phase(this);
	    P2SimpleWrapper P2HardDisk1 = new P2SimpleWrapper(this,new PotentialHardDisk(this));
	    Controller controller1 = new Controller(this);
	    DisplayPhase displayPhase1 = new DisplayPhase(this);
	    IntegratorMD.Timer timer = integratorHard1.new Timer(integratorHard1.chronoMeter());
	    timer.setUpdateInterval(10);
		setBackground(java.awt.Color.yellow);
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        java.awt.Frame f = new java.awt.Frame();   //create a window
        f.setSize(600,350);
        
        Simulation sim = new HSMD2D();
		sim.elementCoordinator.go(); 
		
        f.add(sim.panel());
        
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }//end of main
    
}