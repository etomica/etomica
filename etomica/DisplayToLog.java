package etomica;
import java.io.*;

public class DisplayToLog extends Display {
    
    private PrintWriter logFile;
    
    public DisplayToLog() {
        this(Simulation.instance);
    }
    public DisplayToLog(Simulation sim) {
        this(sim, "default.log");
    }
    public DisplayToLog(String fileName) {
        this(Simulation.instance, fileName);
    }
    public DisplayToLog(Simulation sim, String fileName) {
        super(sim);
        try {
            logFile = new PrintWriter(
                        new BufferedWriter(
                            new FileWriter(fileName)));
        }
        catch(IOException ioe) {
            System.err.println("I/O Exception in DisplayToLog");
    	    System.err.println("Cannot open file: " + fileName);
        }
    }
    
    public java.awt.Component graphic(Object obj) {return null;}
    
    public void intervalAction(Integrator.IntervalEvent evt) {
        if(evt.type() == Integrator.IntervalEvent.START) writeHeading();
        else if(evt.type() == Integrator.IntervalEvent.INTERVAL && --iieCount == 0) {
	        iieCount = updateInterval;
	        doUpdate();
	    }
	    else if(evt.type() == Integrator.IntervalEvent.DONE) writeEnding();
    }
    
    private void writeHeading() {
        logFile.println("Heading");
    }
    
    private void writeEnding() {
        logFile.println("Ending");
        logFile.close();//should do this when controller announces end
    }
    
    public void doUpdate() {
        logFile.println("Update");
    }
    
    //tests and demonstrates use of this class
    public static void main(String[] args) {
        
        Default.ATOM_SIZE = 1.0;
        Simulation sim = new Simulation(new Space2D());
        Simulation.instance = sim;
        
        Species species = new SpeciesDisks(32);
        Potential potential = new PotentialHardDisk();
        Potential2 p2 = new P2SimpleWrapper(potential);
        Integrator integrator = new IntegratorHard();
        Controller controller = new Controller();
        Meter meter = new MeterPressureHard();
        DisplayBox box = new DisplayBox();
        Phase phase = new Phase();
        DisplayPhase display = new DisplayPhase();
        
        DisplayToLog displayLog = new DisplayToLog();
                
        box.setMeter(meter);
        box.setPrecision(8);
        
        sim.elementCoordinator.go();
        phase.setDensity(0.5);
        display.setScale(2.0);
        controller.setMaxSteps(100);

        java.awt.Frame f = new java.awt.Frame();   //create a window
        f.setSize(600,350);
        f.add(sim);         //access the static instance of the simulation to
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }
}