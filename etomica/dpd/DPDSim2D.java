/*
 * Created on Apr 22, 2003
 *
 */
package etomica.dpd;
import etomica.Controller;
import etomica.Default;
import etomica.Phase;
import etomica.Simulation;
import etomica.SpeciesSpheresMono;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterTemperature;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPhase;
import etomica.graphics.SimulationGraphic;
import etomica.nbr.cell.IteratorFactoryCell;
import etomica.space2d.Space2D;




/**
 * Class that contains the main method for running a dissipative particle
 * dynamics (DPD) simulation.
 * 
 * @author cribbin
 * @deprecated
 *
 */

/*Patterned on SwMd2D*/
public class DPDSim2D extends SimulationGraphic{

	public IntegratorDPD integrator;
	public SpeciesSpheresMono species;
	public Phase phase;
	public P2DPD potential;
	public Controller controller;
	public DisplayPhase display;

//	Constructor 1
	public DPDSim2D(){
		/**
		 * @deprecated
		 */
		super(new Space2D());	//creates the simulation and its space. 
		setIteratorFactory(new IteratorFactoryCell(this));
		Simulation.instance = this;
		Default.ATOM_SIZE = 2.0;
		integrator = new IntegratorDPD(this);
		integrator.setInterval(5);
		integrator.setSleepPeriod(1);
		
		integrator.setTimeStep(0.04);	//From Groot & Warren  The constructor does this anyway.
		//integrator.setTemperature(300.);
		
		species = new SpeciesSpheresMono(this);
		species.setNMolecules(400);
		phase = new Phase(this);
		potential = new P2DPD();
		potential.setSpecies(species, species);
		controller = new Controller(this);
		display = new DisplayPhase(this);		//makes the pretty box

		panel().setBackground(java.awt.Color.gray);	//set the background color.
		display.setColorScheme(new etomica.graphics.ColorSchemeRandom());
		
		MeterTime time = new MeterTime(this);
		DisplayBox boxTime = new DisplayBox(this);
		boxTime.setMeter(time);
		boxTime.setPrecision(10);

		MeterEnergy energy = new MeterEnergy(this);
		DisplayBox boxEnergy = new DisplayBox(this);
		boxEnergy.setMeter(energy);
		boxEnergy.setPrecision(10);		
		
		MeterTemperature temperature = new MeterTemperature(this);
		DisplayBox boxTemperature = new DisplayBox(this);
		boxTemperature.setMeter(temperature);
		boxTemperature.setPrecision(10);
				
	}//End-Constructor 1
	

	
	
	
	/**
	 * Runs the class and demonstrates how it is implemented.
	 */
	public static void main(String[] args) {
		DPDSim2D sim = new DPDSim2D();
//		sim.elementCoordinator.go();
		sim.makeAndDisplayFrame();
	}//end method-psvm
}//end class-DPDSim2D
