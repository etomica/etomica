package etomica.modules.dcvgcmd;

import java.awt.Color;
import java.awt.GridBagConstraints;

import javax.swing.JPanel;
import javax.swing.border.TitledBorder;

import etomica.action.Action;
import etomica.action.ActionGroupSeries;
import etomica.action.SimulationRestart;
import etomica.atom.AtomFilter;
import etomica.atom.AtomTypeGroup;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.DataTableAverages;
import etomica.data.meter.MeterNMolecules;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceToggleButton;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTable;
import etomica.graphics.DisplayTextBox;
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
        getDisplayBox(sim.box).setPixelUnit(new Pixel(10));

        getController().getDataStreamPumps().add(sim.profile1pump);
        getController().getDataStreamPumps().add(sim.profile2pump);
        
	    Color colorA = Color.blue;
	    Color colorB = Color.white;

	    GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

	    //Button for cutaway view
	    CutAway cutawayFilter = new CutAway();
	    getDisplayBox(sim.box).setAtomFilter(cutawayFilter);
	    DeviceToggleButton cutawayButton = new DeviceToggleButton(sim.getController());
	    cutawayButton.setModifier(cutawayFilter, "Restore", "Cut tube");
	    cutawayButton.setPostAction(getDisplayBoxPaintAction(sim.box));

	    //Number of each type of atom
	    MeterNMolecules meterA = new MeterNMolecules();
	    MeterNMolecules meterB = new MeterNMolecules();
	    meterA.setBox(sim.box);
	    meterA.setSpecies(sim.species);
	    meterB.setBox(sim.box);
	    meterB.setSpecies(sim.species1);
	    DisplayTextBox boxA = new DisplayTextBox(meterA.getDataInfo());
	    DisplayTextBox boxB = new DisplayTextBox(meterB.getDataInfo());
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
		DisplayTextBox box1 = new DisplayTextBox(sim.thermometer.getDataInfo());
	    final DataPump tpump = new DataPump(sim.thermometer, box1);
        sim.integratorDCV.addIntervalAction(tpump);
		sim.integratorDCV.setActionInterval(tpump, 100);
	    box1.setUnit((Kelvin.UNIT));
		temperatureSlider.setPostAction(new Action() {
			public void actionPerformed() {
				tpump.actionPerformed();
			}
		});

		// Data table tab page
	    DataTableAverages dataTable = new DataTableAverages(sim.integratorDCV.integratormd);
	    dataTable.addDataSource(sim.meterFlux0);
	    dataTable.addDataSource(sim.meterFlux1);
	    dataTable.addDataSource(sim.meterFlux2);
	    dataTable.addDataSource(sim.meterFlux3);
	    DisplayTable table = new DisplayTable(dataTable);
	    
	    table.setTransposed(true);
	    table.setShowingRowLabels(true);
	    table.setRowLabels(new String[] {"Current","Average","Error"});
	    
	    
	    table.setPrecision(7);
	    getPanel().tabbedPane.add("Flux Data", table.graphic());
	    
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
	    ColorSchemeByType colorScheme = (ColorSchemeByType)(getDisplayBox(sim.box).getColorScheme());
		colorScheme.setColor(sim.species.getMoleculeType(),colorA);
		colorScheme.setColor(sim.species1.getMoleculeType(),colorB);
		colorScheme.setColor(((AtomTypeGroup)sim.speciesTube.getMoleculeType()).getChildTypes()[0],java.awt.Color.cyan);

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
	            sim.box.setNMolecules(sim.species, 20);
	            sim.box.setNMolecules(sim.species1, 20);
	            meterAPump.actionPerformed();
	            meterBPump.actionPerformed();
	            tpump.actionPerformed();
	        }
	    });
	    reinitActions.addAction(simRestart);

	    getController().getReinitButton().setAction(reinitActions);
	    getController().getReinitButton().setPostAction(getDisplayBoxPaintAction(sim.box));

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
