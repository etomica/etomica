/*
 * Created on Mar 24, 2005
 */
package etomica.modules.dcvgcmd;

import javax.swing.JPanel;

import etomica.Modifier;
import etomica.action.PhaseImposePbc;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceTrioControllerButton;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayPhaseCanvas3DOpenGL;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTable;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntervalActionAdapter;
import etomica.modules.dcvgcmd.IntegratorDCVGCMD.Mu1Modulator;
import etomica.modules.dcvgcmd.IntegratorDCVGCMD.Mu2Modulator;
import etomica.units.Kelvin;

/**
 * @author msellers and nsives
 *
 */
public class DCVGCMDGraphic extends SimulationGraphic{

	public DCVGCMDGraphic(DCVGCMD sim){
			
	super(sim);	
	
	DisplayPhase display = new DisplayPhase(sim.phase);
    display.setPhase(sim.phase);
    DeviceTrioControllerButton device = new DeviceTrioControllerButton(sim);
    sim.integratorDCV.display = (DisplayPhaseCanvas3DOpenGL)display.graphic();
    
    
    //integratorDCV.setMu(500., 500.);

	DisplayPlot profilePlot = new DisplayPlot();
	add(profilePlot);
	
//Slider to adjust temperature
	DeviceSlider temperatureSlider = new DeviceSlider(sim.getController(), sim.integratorDCV, "temperature");
	temperatureSlider.setUnit((Kelvin.UNIT));
	temperatureSlider.setMinimum(50);
	temperatureSlider.setMaximum(500);
    temperatureSlider.setLabel("Temperature");
			
//Mu Slider Stuff
	Modifier mu1Mod = sim.integratorDCV.new Mu1Modulator(); 
	Modifier mu2Mod = sim.integratorDCV.new Mu2Modulator();
	DeviceSlider mu1Slider = new DeviceSlider(sim.getController(), mu1Mod);
	mu1Slider.setMinimum(-2500);
	mu1Slider.setMaximum(2500);
	DeviceSlider mu2Slider = new DeviceSlider(sim.getController(),mu2Mod);
	mu2Slider.setMinimum(-2500);
	mu2Slider.setMaximum(2500);

//Display to see adjusted temperature
	DisplayBox box1 = new DisplayBox();
    DataPump tpump = new DataPump(sim.thermometer, box1);
	IntervalActionAdapter interval1 = new IntervalActionAdapter (tpump, sim.integratorDCV);
	interval1.setActionInterval(10);
    box1.setUnit((Kelvin.UNIT));
			
	PhaseImposePbc imposePbc = new PhaseImposePbc(sim.phase);
	sim.integratorDCV.addIntervalListener(imposePbc);
	
	DisplayTable table = new DisplayTable();
	add(table);
    sim.fluxAccumulator.addDataSink(table.makeDataSink());
	
	sim.accumulator1.makeDataPusher(new AccumulatorAverage.Type[]{AccumulatorAverage.AVERAGE}).addDataSink(profilePlot.makeDataSink());
	sim.accumulator2.makeDataPusher(new AccumulatorAverage.Type[]{AccumulatorAverage.AVERAGE}).addDataSink(profilePlot.makeDataSink());
	//phase.setDensity(0.5);
//set color of molecules	
	ColorSchemeByType.setColor(sim.species,java.awt.Color.blue);
	ColorSchemeByType.setColor(sim.species1,java.awt.Color.white);
	
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
	 controlPanel.add(temperaturePanel,gbc3);
	 gbc3.gridy = 2;
	 controlPanel.add(muPanel,gbc3);
	
	 panel().remove(panel().devicePanel);
	 
	panel().add(controlPanel);
	 
//panel for membrane choice
	
	
 } //End of constructor

	
	public static void main(String[] arg ){
		
		DCVGCMD sim = new DCVGCMD();
		DCVGCMDGraphic graphic = new DCVGCMDGraphic(sim);
		graphic.makeAndDisplayFrame();
	}//end of main
	
}
