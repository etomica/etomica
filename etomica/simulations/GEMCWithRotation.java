package etomica.simulations;
import etomica.*;
import etomica.units.*;

/**
 * Simple Gibbs-ensemble Monte Carlo simulation of rotating molecules.
 */
//need to update to include setPhaseIteratorFactory
public class GEMCWithRotation extends Simulation {
    
    public GEMCWithRotation() {
        super(new Space2D()/*new Space2DCell()*/);
        Default.ATOM_SIZE = 1.2;
        Default.TEMPERATURE = LennardJones.Temperature.UNIT.toSim(0.420);
        Simulation.instance = this;
        setUnitSystem(new UnitSystem.LJ());
        IntegratorGEMC integratorGEMC1 = new IntegratorGEMC();
	    integratorGEMC1.setDoSleep(false);
	    integratorGEMC1.setInterval(400);
	 //   integratorGEMC1.getMCMoveVolumeExchange().setFrequency(0);
	    
	    SpeciesSpheresRotating speciesDisk1 = new SpeciesSpheresRotating(200);

	    Phase phase1 = new Phase();	    
	    MeterDensity meter1 = new MeterDensity();
	    DisplayBox box1 = new DisplayBox();
	    box1.setMeter(meter1);
	    box1.setUseCurrentValue(false);
	    MCMoveRotate mcRotate1 = new MCMoveRotate(integratorGEMC1);
	    mcRotate1.setPhase(phase1);

	    Phase phase2 = new Phase();
	    MeterDensity meter2 = new MeterDensity();	    
	    DisplayBox box2 = new DisplayBox();
	    box2.setMeter(meter2);
	    box2.setUseCurrentValue(false);
	    MCMoveRotate mcRotate2 = new MCMoveRotate(integratorGEMC1);
	    mcRotate2.setPhase(phase2);
	    
	    P2HardAssociationCone potential = new P2HardAssociationCone();
	    potential.setSigma(speciesDisk1.getDiameter());
	    potential.setWellCutoff(potential.getSigma());
	    Controller controller1 = new Controller();
	    //Configuration displays for each phase
	    DisplayTable table1 = new DisplayTable();
	      
	    DisplayPhase displayPhase1 = new DisplayPhase();
	    DisplayPhase displayPhase2 = new DisplayPhase();
	    
	    //Meters and displays for density in each phase
	    
	    //Slider to adjust temperature
	    DeviceSlider temperatureSlider = new DeviceSlider(integratorGEMC1, "temperature");
	    temperatureSlider.setUnit(new Unit(Kelvin.UNIT));
	    temperatureSlider.setMinimum(50);
	    temperatureSlider.setMaximum(500);

	    integratorGEMC1.setTemperature(Default.TEMPERATURE);
		Simulation.instance.elementCoordinator.go(); 
		
	    /*
	    ColorScheme color1 = new ColorScheme.Simple(java.awt.Color.blue);
	    ColorScheme color2 = new ColorScheme.Simple(java.awt.Color.darkGray);
	    displayPhase1.setColorScheme(color1);
	    displayPhase2.setColorScheme(color2);
	    */
	    phase2.setDensity(0.1);
	    panel().setBackground(java.awt.Color.blue);		
//		((Space2DCell.CellListIteratorFactory)phase1.iteratorFactory()).setNeighborDistance(1.2*Default.ATOM_SIZE);
//        ((Space2DCell.CellListIteratorFactory)phase1.iteratorFactory()).setNCells(10,10);
    }//end of constructor        
        
    public static void main(String[] args) {
        java.awt.Frame f = new java.awt.Frame();   //create a window
        f.setSize(600,350);
        
        Simulation sim = new GEMCWithRotation();
		sim.elementCoordinator.go(); 
        f.add(sim.panel());
        
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
		((Controller)sim.controller(0)).start();
    }//end of main
}