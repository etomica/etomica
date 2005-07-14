package etomica.simulations;

import etomica.Controller;
import etomica.DataSource;
import etomica.Default;
import etomica.Meter;
import etomica.Phase;
import etomica.SpeciesSpheresMono;
import etomica.action.PhaseImposePbc;
import etomica.data.DataSourceAcceptanceRatio;
import etomica.data.DataSourceCountSteps;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceButton;
import etomica.graphics.DeviceTrioControllerButton;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTable;
import etomica.graphics.MediatorGraphic;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorPT;
import etomica.integrator.IntegratorPT.MeterPhaseTracker;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.potential.P2LennardJones;
import etomica.space2d.Space2D;

/**
 * @author David Kofke
 * Simulation demonstrating the use of the ParallelTempering integrator.
 * Simulates a system of Lennard-Jones particles over a range of temperatures.
 * Temperatures are selected such that successive values form a fixed ratio.
 * Keeps averages of acceptance fraction for phase swaps, and presents a plot
 * showing the path followed by each configuration through the phases.
 */
public class ParallelTempering extends SimulationGraphic {

	Phase[] phase;
	SpeciesSpheresMono species;
	Controller controller;
	IntegratorMC[] integrator;
	P2LennardJones potential;
	IntegratorPT.MeterPhaseTracker phaseTracker;


	/**
	 * Constructor for ParallelTempering.
	 * @param nM  number of molecules in each phase
	 * @param nPhase number of phases
	 * @param tMin minimum temperature (LJ units)
	 * @param tMax maximum temperature (LJ units)
	 */
	public ParallelTempering(int nM, int nPhase, double tMin, double tMax) {
		super(Space2D.getInstance());
		Default.makeLJDefaults();
		Default.HISTORY_PERIOD = 1000; //for phase tracking
		         
		elementCoordinator.addMediatorPair(new MediatorGraphic.DisplayMeter.NoAction(elementCoordinator));
        
		species = new SpeciesSpheresMono();
		species.setNMolecules(nM);
		potential = new P2LennardJones();
		potential.setSpecies(species);
        
		IntegratorPT integratorPT = new IntegratorPT(this);
		integratorPT.setDoSleep(true);
		integratorPT.setSleepPeriod(1);
		integratorPT.setInterval(species.getNMolecules());
		integratorPT.setSwapInterval(1000);
        
		phase = new Phase[nPhase];
		integrator = new IntegratorMC[nPhase];
		DataSourceAcceptanceRatio[] meterAccept = new DataSourceAcceptanceRatio[nPhase-1];
		double tLast = 0.0;
		for(int i=0; i<nPhase; i++) {
			double t; 
			t = tMin * Math.exp((double)i/(double)(nPhase-1)*Math.log(tMax/tMin));
            
			phase[i] = new Phase();
			integrator[i] = new IntegratorMC(integratorPT);
			phase[i].setIntegrator(integrator[i]);
			integrator[i].setTemperature(t);
			integrator[i].addListener(new PhaseImposePbc());
            
			MCMoveAtom moveAtom = new MCMoveAtom(integrator[i]);
            
			integratorPT.addIntegrator(integrator[i]);
 
			DisplayPhase display = new DisplayPhase();
			display.setPhase(phase[i]);
			display.setLabel(Double.toString(t));
            
			if(i>0) {
				meterAccept[i-1] = new DataSourceAcceptanceRatio(this, integratorPT.swapMoves()[i-1]);
				meterAccept[i-1].setLabel(Double.toString(tLast)+"<-->"+Double.toString(t));
			}
			tLast = t;
		}
		 
		phaseTracker = integratorPT.new MeterPhaseTracker();
		phaseTracker.setPhases(phase);
		phaseTracker.setUpdateInterval(10);
		DisplayPlot phaseTrackPlot = new DisplayPlot();
		phaseTrackPlot.setDoLegend(false);
		phaseTrackPlot.setWhichValue(MeterAbsMeter);
		phaseTrackPlot.setDataSources(phaseTracker.dataSource());
		phaseTrackPlot.setLabel("Track");
		phaseTrackPlot.setUpdateInterval(10);
        
		DisplayTable acceptanceTable = new DisplayTable();
		acceptanceTable.setDatumSources(meterAccept);
		acceptanceTable.setLabel("Accept");
		acceptanceTable.setWhichValues(new Meter.ValueType[] {Meter.CURRENT, Meter.AVERAGE, Meter.ERROR});
                                
		controller = new Controller();
        
		DisplayBox cyclesBox = new DisplayBox((DatumSource)new DataSourceCountSteps());
		cyclesBox.setPrecision(6);
        
		ColorSchemeByType.setColor(species,java.awt.Color.red);
		DeviceTrioControllerButton buttons = new DeviceTrioControllerButton();
		
//		DeviceButton writeButton = new DeviceButton(this.new WriteHistory());
		WriteHistory writeHistory = this.new WriteHistory();
		writeHistory.setFileName("PT_history");
		writeHistory.setCloseFileEachTime(true);
		DeviceButton writeButton = new DeviceButton(writeHistory.makeWriteAction());
		elementCoordinator.go();
	}

	public static void main(String[] args) {
		int nM = 64;
		int nPhase = 4;
		double tMin = 0.7;
		double tMax = 2.0;

		ParallelTempering ptSim = new ParallelTempering(nM, nPhase, tMin, tMax);
		ptSim.makeAndDisplayFrame();
	}
	
	class WriteHistory extends etomica.log.LoggerAbstract {
//		public WriteHistory() {
//			setLabel("Write");
//		}
		//EXAMPLE
		//	protected void write() throws java.io.IOException {
		//		fileWriter.write(Double.toString(meter.average()));
		//	}
		protected void write() throws java.io.IOException {
			DataSource[] histories = phaseTracker.dataSource();
			int nPhase = histories.length;
			int nPoints = histories[0].data(null).length;
			for(int j=0; j<nPhase; j++) {fileWriter.write(phase[j].integrator().getTemperature()+"\t");}
			fileWriter.write("\n");
			for(int i=0; i<nPoints; i++) {
				if(Double.isNaN(histories[0].data(null)[i])) break;
				for(int j=0; j<nPhase; j++) fileWriter.write((int)histories[j].data(null)[i]+"\t");
				fileWriter.write("\n");
			}
			fileWriter.write("\n");
			
			for(int j=0; j<nPhase; j++) {System.out.print(phase[j].integrator().getTemperature()+"\t");}
			System.out.println();
			for(int i=0; i<nPoints; i++) {
				if(Double.isNaN(histories[0].data(null)[i])) break;
				for(int j=0; j<nPhase; j++) System.out.print((int)histories[j].data(null)[i]+"\t");
				System.out.println();
			}
			System.out.println();
		}
	}
}

