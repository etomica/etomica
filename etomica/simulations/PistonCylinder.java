package etomica.simulations;
import etomica.*;
import etomica.units.*;

/**
 * Simple hard-sphere MD in piston-cylinder apparatus
 */
 
public class PistonCylinder extends Simulation {
    
    public IntegratorHardField integrator;
    public SpeciesSpheres species;
    public SpeciesPistonCylinder speciesPC;
    public Phase phase;
    public P2HardSphere potential;
    public P2HardSphereWall potentialHardDiskWall;
    public Controller controller;
    public DisplayPhase display;
    public MeterCycles meterCycles;
    public DisplayBox displayCycles;

    public PistonCylinder() {
        super(new Space2D());
        Simulation.instance = this;
        
	    species = new SpeciesSpheres(this);
	    speciesPC = new SpeciesPistonCylinder(this);
        speciesPC.setLength(20.);
	    
	    phase = new Phase(this);
	    
	    potential = new P2HardSphere();
	    
	    potentialHardDiskWall = new P2HardSphereWall();
	    
	    integrator = new IntegratorHardField(this);
	    controller = new Controller(this);
	    
	    display = new DisplayPhase(this);
	    
	    meterCycles = new MeterCycles(this);
	    displayCycles = new DisplayBox(this,meterCycles);
	    
		panel().setBackground(java.awt.Color.yellow);

        phase.setBoundary(speciesPC.new Boundary(phase)); //have piston-cylinder system define boundary of phase
        
        //part unique to this class
        etomica.Meter thermometer = new MeterTemperature();
        DisplayBox tBox = new DisplayBox();
        tBox.setMeter(thermometer);
        tBox.setUnit(new Unit(Kelvin.UNIT));
        display.setAlign(1,DisplayPhase.BOTTOM);
        
	    IntegratorMD.Timer timer = integrator.new Timer(integrator.chronoMeter());
	    timer.setUpdateInterval(10);
		Simulation.instance.elementCoordinator.go();

/*        DeviceSlider pressureSlider = new DeviceSlider((SpeciesPistonCylinder.PistonPressureField)phase.firstField(),"pressure");
        Simulation.instance.add(pressureSlider.graphic(null));
        pressureSlider.setUnit(new Unit(Bar.UNIT));
        pressureSlider.setMinimum(50);
        pressureSlider.setMaximum(1000);
*/
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        javax.swing.JFrame f = new javax.swing.JFrame();   //create a window
        f.setSize(600,350);
        
        Simulation sim = new PistonCylinder();
		sim.elementCoordinator.go(); 
		
        f.getContentPane().add(sim.panel());
        
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }//end of main
    
}