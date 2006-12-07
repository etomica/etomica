/*
 * Created on Mar 24, 2005
 */
package etomica.modules.dcvgcmd;

import java.awt.Color;

import javax.swing.JPanel;

import etomica.action.Action;
import etomica.action.ActionGroupSeries;
import etomica.action.SimulationRestart;
import etomica.atom.Atom;
import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomFilter;
import etomica.atom.AtomLeaf;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.DataTableAverages;
import etomica.data.meter.MeterNMolecules;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceToggleButton;
import etomica.graphics.DeviceTrioControllerButton;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTable;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntervalActionAdapter;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierBoolean;
import etomica.units.Kelvin;

/**
 * @author msellers and nsives
 *
 */
public class DCVGCMDGraphic extends SimulationGraphic{

	public DCVGCMDGraphic(final DCVGCMD sim){

	super(sim);	
	
    Color colorA = Color.blue;
    Color colorB = Color.white;
    
    final DisplayPhase displayPhase = getDisplayPhase(sim.phase);
    DeviceTrioControllerButton device = new DeviceTrioControllerButton(sim);
    
    Action repaintAction = new Action() {
        public void actionPerformed() {
            displayPhase.repaint();
        }
        public String getLabel() {
            return null;
        }
    };
    
    //Button for cutaway view
    CutAway cutawayFilter = new CutAway();
    displayPhase.setAtomFilter(cutawayFilter);
    DeviceToggleButton cutawayButton = new DeviceToggleButton(sim.getController());
    cutawayButton.setModifier(cutawayFilter, "Restore", "Cut tube");
    cutawayButton.setPostAction(repaintAction);
    
    
    //integratorDCV.setMu(500., 500.);

	DisplayPlot profilePlot = new DisplayPlot();
    profilePlot.setLabel("Density profile");
    profilePlot.getPlot().setTitle("Density profile");
	add(profilePlot);
    
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
    DataPump meterAPump = new DataPump(meterA,boxA);
    DataPump meterBPump = new DataPump(meterB,boxB);
    new IntervalActionAdapter(meterAPump, sim.integratorDCV);
    new IntervalActionAdapter(meterBPump, sim.integratorDCV);
    JPanel nMoleculePanel = new JPanel();
    nMoleculePanel.add(boxA.graphic());
    nMoleculePanel.add(boxB.graphic());
    nMoleculePanel.setBorder(new javax.swing.border.TitledBorder("Number of atoms"));
	
//Slider to adjust temperature
	DeviceSlider temperatureSlider = new DeviceSlider(sim.getController(), sim.integratorDCV, "temperature");
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
    DataPump tpump = new DataPump(sim.thermometer, box1);
	IntervalActionAdapter interval1 = new IntervalActionAdapter (tpump, sim.integratorDCV);
	interval1.setActionInterval(10);
    box1.setUnit((Kelvin.UNIT));
				
    DataTableAverages dataTable = new DataTableAverages(sim,sim.integratorDCV);
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
	
	sim.accumulator1.addDataSink(profilePlot.getDataSet().makeDataSink(),
            new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
    sim.accumulator2.addDataSink(profilePlot.getDataSet().makeDataSink(),
            new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});

//set color of molecules
    ColorSchemeByType colorScheme = (ColorSchemeByType)displayPhase.getColorScheme();
	colorScheme.setColor(sim.species.getMoleculeType(),colorA);
	colorScheme.setColor(sim.species1.getMoleculeType(),colorB);
	colorScheme.setColor(((AtomFactoryHomo)sim.speciesTube.getFactory()).getChildFactory().getType(),java.awt.Color.cyan);
	
//panel for the start buttons
	  JPanel startPanel = (JPanel)device.graphic();
      
//panel for the temperature control/display
	  JPanel temperaturePanel = new JPanel(new java.awt.GridBagLayout());
	  temperaturePanel.setBorder(new javax.swing.border.TitledBorder("Temperature (K)"));
	  java.awt.GridBagConstraints gbc1 = new java.awt.GridBagConstraints();
	  gbc1.gridx = 0;  gbc1.gridy = 1;
	  gbc1.gridwidth = 1;
	  temperaturePanel.add(temperatureSlider.graphic(null),gbc1);
	  gbc1.gridx = 0;  gbc1.gridy = 0;
	  temperaturePanel.add(box1.graphic(null),gbc1);
	  
//panel for Mu's
	JPanel muPanel = new JPanel(new java.awt.GridBagLayout());
    muPanel.setBorder(new javax.swing.border.TitledBorder("Mu1 and Mu2"));
		java.awt.GridBagConstraints gbc2 = new java.awt.GridBagConstraints();
			gbc2.gridx = 0;  gbc2.gridy = 0;
			gbc2.gridwidth = 1;
			muPanel.add(mu1Slider.graphic(null),gbc2);
			gbc2.gridx = 0;  gbc2.gridy = 1;
			muPanel.add(mu2Slider.graphic(null),gbc2);

	JPanel controlPanel = new JPanel(new java.awt.GridBagLayout());
	 java.awt.GridBagConstraints gbc3 = new java.awt.GridBagConstraints();
	 gbc3.gridy = 0;
	 gbc3.gridx = 0;
	 controlPanel.add(startPanel,gbc3);
	 gbc3.gridy = 1;
     controlPanel.add(cutawayButton.graphic(), gbc3);
     gbc3.gridy = 2;
     controlPanel.add(nMoleculePanel, gbc3);
     gbc3.gridy = 3;
	 controlPanel.add(temperaturePanel,gbc3);
     gbc3.gridy = 4;
	 controlPanel.add(muPanel,gbc3);
	
	 panel().remove(panel().devicePanel);
	 
	panel().add(controlPanel);
    
    SimulationRestart simRestart = (SimulationRestart)device.getReinitButton().getAction();
    simRestart.setConfiguration(sim.config);
    ActionGroupSeries reinitActions = new ActionGroupSeries();
    reinitActions.addAction(new Action() {
        public void actionPerformed() {
            sim.phase.getAgent(sim.species).setNMolecules(20);
            sim.phase.getAgent(sim.species1).setNMolecules(20);
        }
        public String getLabel() {
            return null;
        }
    });
    reinitActions.addAction(simRestart);
    device.getReinitButton().setAction(reinitActions);
    device.getReinitButton().setPostAction(repaintAction);

	 
//	panel for atomsPerRing choice
//	JPanel tubePanel = new JPanel(new java.awt.GridBagLayout());
//	tubePanel.setBorder(new javax.swing.border.TitledBorder("Carbons per Ring"));
//	java.awt.GridBagConstraints gbc4 = new java.awt.GridBagConstraints();
//		gbc4.gridx = 0; gbc4.gridy = 0;
//		gbc4.gridwidth = 1;
//		tubePanel.add(tubeSlider.graphic(null),gbc4);
//		gbc4.gridx = 0;  gbc4.gridy = 1;
//		tubePanel.add(tubeSlider.graphic(null),gbc4);
	
	
 } //End of constructor

	
	public static void main(String[] arg ){
		
		DCVGCMD sim = new DCVGCMD();
        sim.activityIntegrate.setDoSleep(false);
		DCVGCMDGraphic graphic = new DCVGCMDGraphic(sim);
		graphic.makeAndDisplayFrame();
	}//end of main
	
    public static class Applet extends javax.swing.JApplet {

        public void init() {
            getContentPane().add(new DCVGCMDGraphic(new DCVGCMD()).panel());
        }
    }

    private class CutAway implements AtomFilter, ModifierBoolean {
        
        private boolean active = false;
        
        public void setBoolean(boolean b) {active = b;}
        public boolean getBoolean() {return active;}
        
        public boolean accept(Atom atom) {
            if(!active) return true;
            DCVGCMD simulation = (DCVGCMD)getSimulation();
            if(atom.getType().getSpecies() != simulation.speciesTube) return true;
            double x0 = simulation.poreCenter.x(0);
            return ((AtomLeaf)atom).coord.position().x(0) < x0;

        }
    }
}
