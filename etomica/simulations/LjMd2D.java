package etomica.simulations;
import etomica.*;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 2D
 *
 * @author David Kofke
 */
 
public class LjMd2D extends Simulation {
    
    public IntegratorVelocityVerlet integrator;
    public SpeciesDisks species;
    public Phase phase;
    public P2LennardJones potential;
    /*static*/ public Controller controller; //make static for debugging autostart
    public DisplayPhase display;
    public DisplayPlot plot;

    public LjMd2D() {
        super(new Space2D());
        Simulation.instance = this;
	    integrator = new IntegratorVelocityVerlet(this);
	    species = new SpeciesDisks(this);
	    phase = new Phase(this);
	    potential = new P2LennardJones();
	    controller = new Controller(this);
	    display = new DisplayPhase(this);
	    IntegratorMD.Timer timer = integrator.new Timer(integrator.chronoMeter());
	    timer.setUpdateInterval(10);
		setBackground(java.awt.Color.yellow);

		MeterEnergy energy = new MeterEnergy();
		energy.setPhase(phase);
		energy.setHistorying(true);
		energy.setActive(true);		
		energy.getHistory().setNValues(500);		
		plot = new DisplayPlot(this);
		plot.setLabel("Energy");
		plot.setDataSource(energy.getHistory());
		
		integrator.setSleepPeriod(2);
		
		elementCoordinator.go();
        Potential2.Agent potentialAgent = (Potential2.Agent)potential.getAgent(phase);
        potentialAgent.setIterator(new AtomPairIterator(phase));
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {

        Simulation sim = new LjMd2D();
		sim.elementCoordinator.go(); 
		
        javax.swing.JFrame f = new javax.swing.JFrame();   //create a window
        f.setSize(600,350);
        
        f.getContentPane().add(sim.panel());
        
        f.pack();
        f.show();
        f.addWindowListener(Simulation.WINDOW_CLOSER);
     //   controller.start();
    }//end of main
    
}