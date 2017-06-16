/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.dcvgcmd;

import etomica.action.ActionGroupSeries;
import etomica.action.IAction;
import etomica.action.SimulationRestart;
import etomica.atom.AtomFilter;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.data.*;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.meter.MeterNMolecules;
import etomica.graphics.*;
import etomica.listener.IntegratorListenerAction;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierBoolean;
import etomica.molecule.IMolecule;
import etomica.space.Space;
import etomica.units.Kelvin;
import etomica.units.Pixel;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;

/**
 * @author msellers and nsives
 *
 */
public class DCVGCMDGraphic extends SimulationGraphic{

    final static String APP_NAME = "Dual Control-volume GCMD";
    final static int REPAINT_INTERVAL = 70;

    public DCVGCMDGraphic(final DCVGCMD sim, Space _space){

        super(sim, SimulationGraphic.TABBED_PANE, APP_NAME, REPAINT_INTERVAL, _space, sim.getController());	
        getDisplayBox(sim.box).setPixelUnit(new Pixel(7));

        getController().getDataStreamPumps().add(sim.profile1pump);
        getController().getDataStreamPumps().add(sim.profile2pump);
        
        final IAction resetAction = getController().getSimRestart().getDataResetAction();

	    Color[] speciesColors = new Color [] {Color.BLUE, Color.GREEN};

	    GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

	    //Button for cutaway view
	    CutAway cutawayFilter = new CutAway();
	    getDisplayBox(sim.box).setAtomFilter(cutawayFilter);
	    DeviceToggleButton cutawayButton = new DeviceToggleButton(sim.getController());
	    cutawayButton.setModifier(cutawayFilter, "Restore", "Cut tube");
	    cutawayButton.setPostAction(getPaintAction(sim.box));

	    //Number of each type of atom
	    MeterNMolecules meterA = new MeterNMolecules();
	    MeterNMolecules meterB = new MeterNMolecules();
	    meterA.setBox(sim.box);
	    meterA.setSpecies(sim.species1);
	    meterB.setBox(sim.box);
	    meterB.setSpecies(sim.species2);
	    DisplayTextBox boxA = new DisplayTextBox();
	    DisplayTextBox boxB = new DisplayTextBox();
	    boxA.setPrecision(3);
	    boxB.setPrecision(3);
	    boxA.setIntegerDisplay(true);
	    boxB.setIntegerDisplay(true);
	    boxA.setLabel("# Blue Atoms");
	    boxB.setLabel("# Green Atoms");
	    final DataPump meterAPump = new DataPump(meterA,boxA);
	    final DataPump meterBPump = new DataPump(meterB,boxB);
        sim.integratorDCV.getEventManager().addListener(new IntegratorListenerAction(meterAPump));
        sim.integratorDCV.getEventManager().addListener(new IntegratorListenerAction(meterBPump));
	    meterAPump.actionPerformed();
	    meterBPump.actionPerformed();

	    //Slider to adjust temperature
	    DeviceThermoSlider temperatureSlider = new DeviceThermoSlider(sim.getController(), sim.integratorDCV);
		temperatureSlider.setUnit(Kelvin.UNIT);
		temperatureSlider.setMinimum(0);
		temperatureSlider.setMaximum(500);
	    temperatureSlider.setTemperature(Kelvin.UNIT.fromSim(sim.integratorDCV.getTemperature()));
        temperatureSlider.setRadioGroupPostAction(resetAction);

	    //Mu Slider Stuff
		Modifier mu1Mod = sim.integratorDCV.new Mu1Modulator(); 
		Modifier mu2Mod = sim.integratorDCV.new Mu2Modulator();
		DeviceSlider mu1Slider = new DeviceSlider(sim.getController(), mu1Mod);
		mu1Slider.setMinimum(-2500);
		mu1Slider.setMaximum(2500);
		mu1Slider.setShowValues(true);
		mu1Slider.setNMajor(2);
		mu1Slider.setPostAction(resetAction);
		DeviceSlider mu2Slider = new DeviceSlider(sim.getController(),mu2Mod);
		mu2Slider.setMinimum(-2500);
		mu2Slider.setMaximum(2500);
		mu2Slider.setShowValues(true);
		mu2Slider.setNMajor(2);
        mu2Slider.setPostAction(resetAction);

	    //	TubePanel Slider stuff
		//Modifier tubePanelMod = sim.integratorDCV.new tubePanelModifier(); 
		//DeviceSlider tubePanelSlider = new DeviceSlider(sim.getController(), tubePanelMod);
		//tubePanelSlider.setMinimum(8);
		//tubePanelSlider.setMaximum(24);
		
	    //Display to see adjusted temperature
		DisplayTextBox box1 = new DisplayTextBox();
	    final DataPump tpump = new DataPump(sim.thermometer, box1);
	    IntegratorListenerAction tpumpListener = new IntegratorListenerAction(tpump);
        sim.integratorDCV.getEventManager().addListener(tpumpListener);
        tpumpListener.setInterval(100);
	    box1.setUnit((Kelvin.UNIT));
	    box1.setLabel("Measured Temperature");
		temperatureSlider.setSliderPostAction(new IAction() {
			public void actionPerformed() {
				tpump.actionPerformed();
				resetAction.actionPerformed();
			}
		});

		// Data table tab page
        DataTableAverages dataTable = new DataTableAverages(sim.integratorDCV.integratormd, new StatType[]{AccumulatorAverage.AVERAGE, AccumulatorAverage.ERROR},
                1000, new IEtomicaDataSource[0], getController().getDataStreamPumps());
        dataTable.addDataSource(sim.meterFlux0);
        dataTable.addDataSource(sim.meterFlux1);
        dataTable.addDataSource(sim.meterFlux2);
        dataTable.addDataSource(sim.meterFlux3);
        DisplayTable table = new DisplayTable(dataTable);

	    table.setTransposed(true);
	    table.setShowingRowLabels(true);
	    table.setRowLabels(new String[] {"Average","Error"});
	    table.setColumnHeader(new DataTag[]{sim.meterFlux0.getTag()}, "Left (1)");
        table.setColumnHeader(new DataTag[]{sim.meterFlux1.getTag()}, "Right (1)");
        table.setColumnHeader(new DataTag[]{sim.meterFlux2.getTag()}, "Left (2)");
        table.setColumnHeader(new DataTag[]{sim.meterFlux3.getTag()}, "Right (2)");

	    table.setPrecision(7);
	    getPanel().tabbedPane.add("Flux Data", table.graphic());
	    
	    // Density profile tab page
		DisplayPlot profilePlot = new DisplayPlot();
	    profilePlot.setLabel("Density Profile");
	    profilePlot.getPlot().setTitle("Density Profile");
	    profilePlot.getPlot().setColors(speciesColors);
	    profilePlot.setLegend(new DataTag[]{sim.accumulator1.getTag()}, "Density (1)");
        profilePlot.setLegend(new DataTag[]{sim.accumulator2.getTag()}, "Density (2)");
		getPanel().tabbedPane.add("Density Profile", profilePlot.graphic());

		sim.accumulator1.addDataSink(profilePlot.getDataSet().makeDataSink(),
	            new AccumulatorAverage.StatType[]{AccumulatorAverage.AVERAGE});
	    sim.accumulator2.addDataSink(profilePlot.getDataSet().makeDataSink(),
	            new AccumulatorAverage.StatType[]{AccumulatorAverage.AVERAGE});

	    //set color of molecules
	    ColorSchemeByType colorScheme = (ColorSchemeByType)(getDisplayBox(sim.box).getColorScheme());
		colorScheme.setColor(sim.species1.getLeafType(), speciesColors[0]);
		colorScheme.setColor(sim.species2.getLeafType(), speciesColors[1]);
		colorScheme.setColor(sim.speciesTube.getAtomType(0),java.awt.Color.cyan);

		DiameterHashByType diameterHash = (DiameterHashByType)getDisplayBox(sim.box).getDiameterHash();
		diameterHash.setDiameter(sim.species1.getLeafType(), 3.0);
        diameterHash.setDiameter(sim.species2.getLeafType(), 3.0);
        diameterHash.setDiameter(sim.speciesTube.getLeafType(), 3.0);

	    //panel for Mu's
		JPanel muPanel = new JPanel(new java.awt.GridBagLayout());
	    muPanel.setBorder(new TitledBorder(null, "Mu1 and Mu2", TitledBorder.CENTER, TitledBorder.TOP));
		muPanel.add(mu1Slider.graphic(null),vertGBC);
		muPanel.add(mu2Slider.graphic(null),vertGBC);
	
		add(getController());
		add(cutawayButton);
		add(boxA);
		add(boxB);
		add(temperatureSlider);
	    getPanel().controlPanel.add(muPanel,vertGBC);
	    //panel for the temperature control/display
		add(box1);

		
	    SimulationRestart simRestart = (SimulationRestart)getController().getReinitButton().getAction();
	    simRestart.setConfiguration(sim.config);
	    ActionGroupSeries reinitActions = new ActionGroupSeries();
	    reinitActions.addAction(new IAction() {
	        public void actionPerformed() {
	            sim.box.setNMolecules(sim.species1, 20);
	            sim.box.setNMolecules(sim.species2, 20);
	            meterAPump.actionPerformed();
	            meterBPump.actionPerformed();
	            tpump.actionPerformed();
	        }
	    });
	    reinitActions.addAction(simRestart);

	    getController().getReinitButton().setAction(reinitActions);
	    getController().getReinitButton().setPostAction(getPaintAction(sim.box));

		getPanel().toolbar.addContributor("Colin Tedlock");

	} //End of constructor

	
	public static void main(String[] arg ){
		
		DCVGCMD sim = new DCVGCMD();
		DCVGCMDGraphic graphic = new DCVGCMDGraphic(sim, sim.getSpace());
		graphic.makeAndDisplayFrame(APP_NAME);
	}//end of main
	
    public static class Applet extends javax.swing.JApplet {

        public void init() {
    		DCVGCMD sim = new DCVGCMD();
    		DCVGCMDGraphic graphic = new DCVGCMDGraphic(sim, sim.getSpace());
            getContentPane().add(new DCVGCMDGraphic(sim, sim.getSpace()).getPanel());
        }
    }

    private class CutAway implements AtomFilter, ModifierBoolean {
        
        private boolean active = false;
        
        public void setBoolean(boolean b) {active = b;}
        public boolean getBoolean() {return active;}
        
        public boolean accept(IAtom atom) {
            if(!active) return true;
            if(atom.getType().getSpecies() != ((DCVGCMD)simulation).speciesTube) return true;
            double x0 = ((DCVGCMD)simulation).poreCenter.getX(0);
            return atom.getPosition().getX(0) < x0;

        }
        public boolean accept(IMolecule mole) {
            return false;
        }
    }
}
