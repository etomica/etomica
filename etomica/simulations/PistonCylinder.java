package etomica.simulations;
import etomica.*;
import etomica.graphics.*;
import etomica.units.*;

/**
 * Simple hard-sphere MD in piston-cylinder apparatus
 */
public class PistonCylinder extends SimulationGraphic {
    
    public IntegratorHardField integrator;
    public SpeciesSpheresMono species;
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
        
	    species = new SpeciesSpheresMono(this);
	    speciesPC = new SpeciesPistonCylinder(this);
        speciesPC.setLength(20.);
	    
	    phase = new Phase(this);
	    
	    potential = new P2HardSphere();
	    potential.setSpecies(species, species);
	    
	    Potential2Group potentialDiskPC = new Potential2Group();
	    potentialDiskPC.setSpecies(species, speciesPC);
	    
	    potentialHardDiskWall = new P2HardSphereWall(potentialDiskPC, Default.ATOM_SIZE);
	    potentialHardDiskWall.setIterator(new AtomPairIteratorGeneral(Simulation.instance.space,
	                new AtomIteratorSinglet(), new AtomIteratorSequential()));
	    
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
        
	    DisplayTimer timer = new DisplayTimer(integrator);
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
        SimulationGraphic sim = new PistonCylinder();
		sim.elementCoordinator.go();
		sim.makeAndDisplayFrame();
   }//end of main
    
}