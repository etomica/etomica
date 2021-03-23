/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.dcvgcmd;

import etomica.action.ActionGroupSeries;
import etomica.action.IAction;
import etomica.action.SimulationRestart;
import etomica.atom.DiameterHashByType;
import etomica.data.*;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.meter.MeterNMolecules;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.space.Space;
import etomica.units.Kelvin;
import etomica.units.Pixel;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;

/**
 * @author msellers and nsives
 */
public class DCVGCMDGraphicFasterer extends SimulationGraphic {

	final static String APP_NAME = "Dual Control-volume GCMD";
	final static int REPAINT_INTERVAL = 70;

	public DCVGCMDGraphicFasterer(final DCVGCMDFasterer sim, Space _space) {

		super(sim, SimulationGraphic.TABBED_PANE, APP_NAME, REPAINT_INTERVAL);
		getDisplayBox(sim.box).setPixelUnit(new Pixel(7));

		getController().getDataStreamPumps().add(sim.profile1pump);
		getController().getDataStreamPumps().add(sim.profile2pump);

		final IAction resetAction = getController().getSimRestart().getDataResetAction();

		Color[] speciesColors = new Color[]{Color.BLUE, Color.GREEN};

		GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

		//Number of each type of atom
		MeterNMolecules meterA = new MeterNMolecules();
		MeterNMolecules meterB = new MeterNMolecules();
		meterA.setBox(sim.box);
		meterA.setSpecies(sim.propane);
		meterB.setBox(sim.box);
		meterB.setSpecies(sim.propene);
		DisplayTextBox boxA = new DisplayTextBox();
		DisplayTextBox boxB = new DisplayTextBox();
		boxA.setPrecision(3);
		boxB.setPrecision(3);
		boxA.setIntegerDisplay(true);
		boxB.setIntegerDisplay(true);
		boxA.setLabel("# Blue Atoms");
		boxB.setLabel("# Green Atoms");
		final DataPump meterAPump = new DataPump(meterA, boxA);
		final DataPump meterBPump = new DataPump(meterB, boxB);
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
		double muMin = -10000;
		double muMax = 0;
		Modifier muAMMod = new ModifierGeneral(sim.integratorDCV.mcMoves()[0], "mu");
		Modifier muAPMod = new ModifierGeneral(sim.integratorDCV.mcMoves()[1], "mu");
		Modifier muBMMod = new ModifierGeneral(sim.integratorDCV.mcMoves()[2], "mu");
		Modifier muBPMod = new ModifierGeneral(sim.integratorDCV.mcMoves()[3], "mu");
		DeviceSlider muAMSlider = new DeviceSlider(sim.getController(), muAMMod);
		muAMSlider.setMinimum(muMin);
		muAMSlider.setMaximum(muMax);
		muAMSlider.setShowValues(true);
		muAMSlider.setNMajor(2);
		muAMSlider.setPostAction(resetAction);
		muAMSlider.setUnit(Kelvin.UNIT);
		DeviceSlider muAPSlider = new DeviceSlider(sim.getController(), muAPMod);
		muAPSlider.setMinimum(muMin);
		muAPSlider.setMaximum(muMax);
		muAPSlider.setShowValues(true);
		muAPSlider.setNMajor(2);
		muAPSlider.setPostAction(resetAction);
		muAPSlider.setUnit(Kelvin.UNIT);
		DeviceSlider muBMSlider = new DeviceSlider(sim.getController(), muBMMod);
		muBMSlider.setMinimum(muMin);
		muBMSlider.setMaximum(muMax);
		muBMSlider.setShowValues(true);
		muBMSlider.setNMajor(2);
		muBMSlider.setPostAction(resetAction);
		muBMSlider.setUnit(Kelvin.UNIT);
		DeviceSlider muBPSlider = new DeviceSlider(sim.getController(), muBPMod);
		muBPSlider.setMinimum(muMin);
		muBPSlider.setMaximum(muMax);
		muBPSlider.setShowValues(true);
		muBPSlider.setNMajor(2);
		muBPSlider.setPostAction(resetAction);
		muBPSlider.setUnit(Kelvin.UNIT);

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
				1000, new IDataSource[0], getController().getDataStreamPumps());
		dataTable.addDataSource(sim.meterFlux0);
		dataTable.addDataSource(sim.meterFlux1);
		dataTable.addDataSource(sim.meterFlux2);
		dataTable.addDataSource(sim.meterFlux3);
		DisplayTable table = new DisplayTable(dataTable);

		table.setTransposed(true);
		table.setShowingRowLabels(true);
		table.setRowLabels(new String[]{"Average", "Error"});
		table.setColumnHeader(new DataTag[]{sim.meterFlux0.getTag()}, "Left (1)");
		table.setColumnHeader(new DataTag[]{sim.meterFlux1.getTag()}, "Right (1)");
		table.setColumnHeader(new DataTag[]{sim.meterFlux2.getTag()}, "Left (2)");
		table.setColumnHeader(new DataTag[]{sim.meterFlux3.getTag()}, "Right (2)");

		table.setPrecision(7);
		getPanel().tabbedPane.add("Flux Data", table.graphic());

		// Density profile tab page
		DisplayPlotXChart profilePlot = new DisplayPlotXChart();
		profilePlot.setLabel("Density Profile");
		profilePlot.getPlot().setTitle("Density Profile");
		profilePlot.setLegend(new DataTag[]{sim.accumulator1.getTag()}, "Density (1)");
		profilePlot.setLegend(new DataTag[]{sim.accumulator2.getTag()}, "Density (2)");
		getPanel().tabbedPane.add("Density Profile", profilePlot.graphic());

		sim.accumulator1.addDataSink(profilePlot.getDataSet().makeDataSink(),
				new StatType[]{sim.accumulator1.AVERAGE});
		sim.accumulator2.addDataSink(profilePlot.getDataSet().makeDataSink(),
				new StatType[]{sim.accumulator2.AVERAGE});

		profilePlot.getPlot().setColors(speciesColors);

		DisplayPlotXChart temperaturePlot = new DisplayPlotXChart();
		temperaturePlot.setUnit(Kelvin.UNIT);
		temperaturePlot.setLabel("Temperature Profile");
		temperaturePlot.getPlot().setTitle("Temperature Profile");
		temperaturePlot.setLegend(new DataTag[]{sim.profileTemperature1.getTag()}, "Temperature (1)");
		temperaturePlot.setLegend(new DataTag[]{sim.profileTemperature2.getTag()}, "Temperature (2)");
		getPanel().tabbedPane.add("Temperature Profile", temperaturePlot.graphic());

		sim.profile1TemperaturePump.setDataSink(temperaturePlot.getDataSet().makeDataSink());
		sim.profile2TemperaturePump.setDataSink(temperaturePlot.getDataSet().makeDataSink());

		temperaturePlot.getPlot().setColors(speciesColors);

		//set color of molecules
		ColorSchemeByType colorScheme = (ColorSchemeByType) (getDisplayBox(sim.box).getColorScheme());
		colorScheme.setColor(sim.propane.getTypeByName("propaneCH3"), speciesColors[0]);
		colorScheme.setColor(sim.propane.getTypeByName("propaneCH2"), speciesColors[0]);
		colorScheme.setColor(sim.propene.getTypeByName("propeneCH3"), speciesColors[1]);
		colorScheme.setColor(sim.propene.getTypeByName("propeneCH2"), speciesColors[1]);
		colorScheme.setColor(sim.propene.getTypeByName("propeneCH"), speciesColors[1]);
		colorScheme.setColor(sim.membrane.getAtomType(0), Color.cyan);

		DiameterHashByType diameterHash = (DiameterHashByType) getDisplayBox(sim.box).getDiameterHash();
		diameterHash.setDiameter(sim.propane.getTypeByName("propaneCH3"), 3.75);
		diameterHash.setDiameter(sim.propane.getTypeByName("propaneCH2"), 3.95);
		diameterHash.setDiameter(sim.propene.getTypeByName("propeneCH3"), 3.75);
		diameterHash.setDiameter(sim.propene.getTypeByName("propeneCH2"), 3.675);
		diameterHash.setDiameter(sim.propene.getTypeByName("propeneCH"), 3.73);
		diameterHash.setDiameter(sim.membrane.getAtomType(0), sim.sigma);

		//panel for Mu's
		JPanel muPanel = new JPanel(new GridBagLayout());
		muPanel.setBorder(new TitledBorder(null, "Mu1 and Mu2 (K)", TitledBorder.CENTER, TitledBorder.TOP));
		muPanel.add(muAMSlider.graphic(), vertGBC);
		muPanel.add(muAPSlider.graphic(), vertGBC);
		muPanel.add(muBMSlider.graphic(), vertGBC);
		muPanel.add(muBPSlider.graphic(), vertGBC);

		add(getController());
		add(boxA);
		add(boxB);
		add(temperatureSlider);
		getPanel().controlPanel.add(muPanel, vertGBC);
		//panel for the temperature control/display
		add(box1);


		SimulationRestart simRestart = (SimulationRestart) getController().getReinitButton().getAction();
		ActionGroupSeries reinitActions = new ActionGroupSeries();
		reinitActions.addAction(new IAction() {
			public void actionPerformed() {
				sim.box.setNMolecules(sim.propane, 20);
				sim.box.setNMolecules(sim.propene, 20);
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


	public static void main(String[] arg) {

		DCVGCMDFasterer sim = new DCVGCMDFasterer();
		DCVGCMDGraphicFasterer graphic = new DCVGCMDGraphicFasterer(sim, sim.getSpace());
		graphic.makeAndDisplayFrame(APP_NAME);
	}//end of main

}
