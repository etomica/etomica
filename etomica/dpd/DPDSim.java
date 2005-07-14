/*
 * Created on Apr 22, 2003
 *
 */
package etomica.dpd;
import etomica.Controller;
import etomica.Default;
import etomica.Phase;
import etomica.Simulation;
import etomica.Space;
import etomica.SpeciesSpheresMono;
import etomica.data.DataSourceCountTime;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterTemperature;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceTrioControllerButton;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayPhaseCanvas3DOpenGL;
import etomica.graphics.SimulationGraphic;
import etomica.nbr.cell.IteratorFactoryCell;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;


/**
 * Class that contains the main method for running a dissipative particle
 * dynamics (DPD) simulation.
 * 
 * @author cribbin
 * 
 *  
 */

/*Patterned on SwMd2D*/
public class DPDSim extends SimulationGraphic{

	public IntegratorDPD integrator;
	public SpeciesSpheresMono species;
	public Phase phase;
	public P2DPD potential;
	public Controller controller;
	public DisplayPhase display;

/**
 * Default constructor, making a 3-d simulation
 * @see java.lang.Object#Object()
 */
	public DPDSim() {this(3);}
	
	public DPDSim(int dim){
		super(dim==2 ? (Space)Space2D.getInstance(): (Space)Space3D.getInstance());	//creates the simulation and its space. 
		Default.makeLJDefaults();
		double density = 3.0;
		double temperature = Default.TEMPERATURE;
		setIteratorFactory(new IteratorFactoryCell(this,6));
		Simulation.instance = this;
		integrator = new IntegratorDPD(this);
		integrator.setInterval(1);
		integrator.setSleepPeriod(5);
		
		integrator.setTimeStep(0.04);	//From Groot & Warren  The constructor does this anyway.
		integrator.setLambdaV(0.65);
		species = new SpeciesSpheresMono(this);
		species.setNMolecules(512);
		phase = new Phase(this);
		potential = new P2DPD(3.0,75.0/density,1.0);
//		potential = new P2DPD(3.0,1.0,1.0);
		potential.setTemperature(temperature);
		potential.setSpecies(species, species);
		controller = new Controller(this);
		DeviceTrioControllerButton controllerButton = new DeviceTrioControllerButton(this);

		panel().setBackground(java.awt.Color.gray);	//set the background color.
//		display.setColorScheme(new etomica.graphics.ColorSchemeRandom());
		ColorSchemeByType.setColor(species, java.awt.Color.red);

		DataSourceCountTime time = new DataSourceCountTime(integrator);
		DisplayBox boxTime = new DisplayBox(this);
		boxTime.setMeter(time);
		boxTime.setPrecision(10);

		MeterKineticEnergy energy = new MeterKineticEnergy(this);
		DisplayBox boxEnergy = new DisplayBox(this);
		boxEnergy.setMeter(energy);
		boxEnergy.setPrecision(10);		
		
		MeterTemperature meterT = new MeterTemperature(this);
		DisplayBox boxTemperature = new DisplayBox(this);
		boxTemperature.setMeter(meterT);
		boxTemperature.setPrecision(10);
				
//		elementCoordinator.go();
		
		DeviceSlider densitySlider = new DeviceSlider(phase, "density");
		densitySlider.setPrecision(2);
		densitySlider.setMaximum(5.0);
		densitySlider.setMinimum(0.0);
		densitySlider.setNMajor(5);
		densitySlider.setValue(density);

		display = new DisplayPhase(this);		//makes the pretty box
//		display.setUpdateInterval(1);
		((DisplayPhaseCanvas3DOpenGL)display.canvas).setDrawExpansionFactor(1.2);
//		elementCoordinator.go();						
	}//End-Constructor 1
	
	/**
	 * Runs the class and demonstrates how it is implemented.
	 */
	public static void main(String[] args) {
		int dim = 3;	//Dimension of simulation  
		DPDSim sim = new DPDSim(dim);
//		sim.elementCoordinator.go();
		sim.makeAndDisplayFrame();
	}//end method-psvm
}//end class-DPDSim3D
