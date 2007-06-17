package etomica.modules.dcvgcmd;

import java.awt.Color;
import java.awt.GridBagConstraints;

import javax.swing.JPanel;
import javax.swing.border.TitledBorder;

import etomica.action.Action;
import etomica.action.ActionGroupSeries;
import etomica.action.SimulationRestart;
import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomFilter;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.DataTableAverages;
import etomica.data.meter.MeterNMolecules;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceToggleButton;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTable;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierBoolean;
import etomica.units.Kelvin;
import etomica.units.Pixel;

/**
 * @author msellers and nsives
 *
 */
public class DCVGCMDGraphic extends SimulationGraphic{

	final static String APP_NAME = "Dual Control-volume GCMD";
	final static int REPAINT_INTERVAL = 70;

	public DCVGCMDGraphic(final DCVGCMD sim){

		super(sim, SimulationGraphic.TABBED_PANE, APP_NAME, REPAINT_INTERVAL);	
        getDisplayPhase(sim.phase).setPixelUnit(new Pixel(10));

	    Color colorA = Color.blue;
	    Color colorB = Color.white;

	    GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

	    //Button for cutaway view
	    CutAway cutawayFilter = new CutAway();
	    getDisplayPhase(sim.phase).setAtomFilter(cutawayFilter);
	    DeviceToggleButton cutawayButton = new DeviceToggleButton(sim.getController());
	    cutawayButton.setModifier(cutawayFilter, "Restore", "Cut tube");
	    cutawayButton.setPostAction(getDisplayPhasePaintAction(sim.phase));

	    //Number of each type of atom
	    MeterNMolecules meterA = new MeterNMolecules();
	    MeterNMolecules meterB = new MeterNMolecules();
	    meterA.setPhase(sim.phase);
	    meterA.setSpecies(sim.species);
	    meterB.setPhase(sim.phase);
	    meterB.setSpecies(sim.species1);
	    DisplayBox boxA = new DisplayBox(meterA.getDataInfo());
	    DisplayBox boxB = new DisplayBox(meterB.getDataInfo());
	    boxA.setPrecision(3);
	    boxB.setPrecision(3);
	    boxA.setIntegerDisplay(true);
	    boxB.setIntegerDisplay(true);
	    boxA.setLabel("  Blue  ");
	    boxB.setLabel(" White  ");
	    final DataPump meterAPump = new DataPump(meterA,boxA);
	    final DataPump meterBPump = new DataPump(meterB,boxB);
        sim.integratorDCV.addIntervalAction(meterAPump);
        sim.integratorDCV.addIntervalAction(meterBPump);
	    JPanel nMoleculePanel = new JPanel();
	    nMoleculePanel.add(boxA.graphic());
	    nMoleculePanel.add(boxB.graphic());
	    nMoleculePanel.setBorder(new TitledBorder("Number of atoms"));
	    meterAPump.actionPerformed();
	    meterBPump.actionPerformed();

	    //Slider to adjust temperature
		final DeviceSlider temperatureSlider = new DeviceSlider(sim.getController(), sim.integratorDCV, "temperature");
		temperatureSlider.setUnit(Kelvin.UNIT);
		temperatureSlider.setMinimum(50);
		temperatureSlider.setMaximum(500);
	    temperatureSlider.setLabel("Temperature");
	    temperatureSlider.setValue(Kelvin.UNIT.fromSim(sim.integratorDCV.getTemperature()));
	    
	    //Mu Slider Stuff
		Modifier mu1Mod = sim.integratorDCV.new Mu1Modulator(); 
		Modifier mu2Mod = sim.integratorDCV.new Mu2Modulator();
		DeviceSlider mu1Slider = new DeviceSlider(sim.getController(), mu1Mod);
		mu1Slider.setMinimum(-2500);
		mu1Slider.setMaximum(2500);
		DeviceSlider mu2Slider = new DeviceSlider(sim.getController(),mu2Mod);
		mu2Slider.setMinimum(-2500);
		mu2Slider.setMaximum(2500);
	
	    //	TubePanel Slider stuff
		//Modifier tubePanelMod = sim.integratorDCV.new tubePanelModifier(); 
		//DeviceSlider tubePanelSlider = new DeviceSlider(sim.getController(), tubePanelMod);
		//tubePanelSlider.setMinimum(8);
		//tubePanelSlider.setMaximum(24);
		
	    //Display to see adjusted temperature
		DisplayBox box1 = new DisplayBox(sim.thermometer.getDataInfo());
	    final DataPump tpump = new DataPump(sim.thermometer, box1);
        sim.integratorDCV.addIntervalAction(tpump);
		sim.integratorDCV.setActionInterval(tpump, 100);
	    box1.setUnit((Kelvin.UNIT));
		temperatureSlider.setPostAction(new Action() {
			public void actionPerformed() {
				tpump.actionPerformed();
			}
		});

	    DataTableAverages dataTable = new DataTableAverages(sim.integratorDCV);
	    dataTable.addDataSource(sim.meterFlux0);
	    dataTable.addDataSource(sim.meterFlux1);
	    dataTable.addDataSource(sim.meterFlux2);
	    dataTable.addDataSource(sim.meterFlux3);
	    DisplayTable table = new DisplayTable(dataTable);
		add(table);
	    table.setRowLabels(new String[] { "Current", "Average", "Error" });
	    table.setTransposed(true);
	    table.setShowingRowLabels(true);
	    table.setPrecision(7);

	    // Density profile tab page
		DisplayPlot profilePlot = new DisplayPlot();
	    profilePlot.setLabel("Density profile");
	    profilePlot.getPlot().setTitle("Density profile");
		getPanel().tabbedPane.add("Density profile", profilePlot.graphic());

		sim.accumulator1.addDataSink(profilePlot.getDataSet().makeDataSink(),
	            new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
	    sim.accumulator2.addDataSink(profilePlot.getDataSet().makeDataSink(),
	            new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});

	    //set color of molecules
	    ColorSchemeByType colorScheme = (ColorSchemeByType)(getDisplayPhase(sim.phase).getColorScheme());
		colorScheme.setColor(sim.species.getMoleculeType(),colorA);
		colorScheme.setColor(sim.species1.getMoleculeType(),colorB);
		colorScheme.setColor(((AtomFactoryHomo)sim.speciesTube.getFactory()).getChildFactory().getType(),java.awt.Color.cyan);

	    //panel for the temperature control/display
		JPanel temperaturePanel = new JPanel(new java.awt.GridBagLayout());
		temperaturePanel.setBorder(new TitledBorder("Temperature (K)"));
		temperaturePanel.add(box1.graphic(null),vertGBC);
		temperaturePanel.add(temperatureSlider.graphic(null),vertGBC);

	    //panel for Mu's
		JPanel muPanel = new JPanel(new java.awt.GridBagLayout());
	    muPanel.setBorder(new TitledBorder("Mu1 and Mu2"));
		muPanel.add(mu1Slider.graphic(null),vertGBC);
		muPanel.add(mu2Slider.graphic(null),vertGBC);
	
		add(getController());
		add(cutawayButton);
	    getPanel().controlPanel.add(nMoleculePanel, vertGBC);
	    getPanel().controlPanel.add(temperaturePanel,vertGBC);
	    getPanel().controlPanel.add(muPanel,vertGBC);
		
	    SimulationRestart simRestart = (SimulationRestart)getController().getReinitButton().getAction();
	    simRestart.setConfiguration(sim.config);
	    ActionGroupSeries reinitActions = new ActionGroupSeries();
	    reinitActions.addAction(new Action() {
	        public void actionPerformed() {
	            sim.phase.getAgent(sim.species).setNMolecules(20);
	            sim.phase.getAgent(sim.species1).setNMolecules(20);
	            meterAPump.actionPerformed();
	            meterBPump.actionPerformed();
	            tpump.actionPerformed();
	        }
	    });
	    reinitActions.addAction(simRestart);

	    getController().getReinitButton().setAction(reinitActions);
	    getController().getReinitButton().setPostAction(getDisplayPhasePaintAction(sim.phase));

		getPanel().toolbar.addContributor("Colin Tedlock");

	} //End of constructor

	
	public static void main(String[] arg ){
		
		DCVGCMD sim = new DCVGCMD();
        sim.activityIntegrate.setDoSleep(false);
		DCVGCMDGraphic graphic = new DCVGCMDGraphic(sim);
		graphic.makeAndDisplayFrame(APP_NAME);
	}//end of main
	
    public static class Applet extends javax.swing.JApplet {

        public void init() {
            getContentPane().add(new DCVGCMDGraphic(new DCVGCMD()).getPanel());
        }
    }

    private class CutAway implements AtomFilter, ModifierBoolean {
        
        private boolean active = false;
        
        public void setBoolean(boolean b) {active = b;}
        public boolean getBoolean() {return active;}
        
        public boolean accept(IAtom atom) {
            if(!active) return true;
            if(atom.getType().getSpecies() != ((DCVGCMD)simulation).speciesTube) return true;
            double x0 = ((DCVGCMD)simulation).poreCenter.x(0);
            return ((IAtomPositioned)atom).getPosition().x(0) < x0;

        }
    }
}
