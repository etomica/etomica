/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.reactionequilibrium;

import etomica.action.ActionGroupSeries;
import etomica.action.IAction;
import etomica.action.SimulationRestart;
import etomica.action.controller.Controller;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.*;
import etomica.data.types.DataTable;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.*;
import etomica.graphics.DisplayTextBox.LabelType;
import etomica.integrator.IntegratorListenerAction;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.modifier.ModifierSQWEpsilon;
import etomica.potential.P2HardGeneric;
import etomica.space.Space;
import etomica.species.SpeciesGeneral;
import etomica.units.Angstrom;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;
import etomica.util.Constants.CompassDirection;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;

/**
 * @author William Scharmach
 */
public class ReactionEquilibriumGraphic extends SimulationGraphic {

	private static final String APP_NAME = "Reaction Equilibrium";
	private static final int REPAINT_INTERVAL = 1;
	protected final ReactionEquilibrium sim;
	protected boolean initializing;
	protected AccumulatorAverage densityAccum;
	protected DataPump dimerPump;
	protected DisplayTextBoxesCAE densityDisplay;
	private DeviceThermoSlider temperatureSelect;
	protected final IAction resetAction;

	public ReactionEquilibriumGraphic(ReactionEquilibrium simulation, Space space) {

		super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);
		this.sim = simulation;

		resetAction = getController().getSimRestart().getDataResetAction();

		GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

		getDisplayBox(sim.box).setPixelUnit(new Pixel(5));

		temperatureSelect = new DeviceThermoSlider(sim.getController(), sim.integratorHard1);
		temperatureSelect.setUnit(Kelvin.UNIT);
		temperatureSelect.setMaximum(2500);
		temperatureSelect.setTemperature(300); //sets 300K as selected temperature
		temperatureSelect.setIsothermal();
		temperatureSelect.setSliderPostAction(resetAction);
		temperatureSelect.setRadioGroupPostAction(resetAction);
		((ColorSchemeByType) getDisplayBox(sim.box).getColorScheme()).setColor(sim.speciesA.getLeafType(), Color.RED);
		((ColorSchemeByType) getDisplayBox(sim.box).getColorScheme()).setColor(sim.speciesB.getLeafType(), Color.BLACK);

		//	adjustment of species properties
		MySpeciesEditor AEditor = new MySpeciesEditor(sim.speciesA, "Red");
		MySpeciesEditor BEditor = new MySpeciesEditor(sim.speciesB, "Black");
		int ms = 10;
		AEditor.nSlider.getSlider().setMajorTickSpacing(ms);
		BEditor.nSlider.getSlider().setMajorTickSpacing(ms);
		AEditor.nSlider.getSlider().setMinorTickSpacing(1);
		BEditor.nSlider.getSlider().setMinorTickSpacing(1);
		AEditor.nSlider.getSlider().setLabelTable(
				AEditor.nSlider.getSlider().createStandardLabels(ms));
		BEditor.nSlider.getSlider().setLabelTable(
				AEditor.nSlider.getSlider().createStandardLabels(ms));

		MassEditor AMassEditor = new MassEditor(sim.getController(), sim.speciesA, "Red");
		MassEditor BMassEditor = new MassEditor(sim.getController(), sim.speciesB, "Black");

		ActionGroupSeries fullResetAction = new ActionGroupSeries(resetAction, () -> sim.integratorHard1.reset());

		//sliders to adjust potentials well depth
		int eMin = 0;
		int eMax = 1000;
		//    final DeviceSlider AASlider = new DeviceSlider(AAbonded, "epsilon");
		final DeviceSlider AASlider = new DeviceSlider(sim.getController(), new ModifierSQWEpsilon(sim.AAbonded));
		AASlider.setUnit((Kelvin.UNIT));
		AASlider.setShowBorder(true);
		AASlider.setLabel("RR epsilon");
		AASlider.setMinimum(eMin);
		AASlider.setMaximum(eMax);
		AASlider.setNMajor(5);
		AASlider.getSlider().setSnapToTicks(true);
		AASlider.setPostAction(fullResetAction);
		//    final DeviceSlider ABSlider = new DeviceSlider(ABbonded, "epsilon");
		final DeviceSlider ABSlider = new DeviceSlider(sim.getController(), new ModifierSQWEpsilon(sim.ABbonded));
		ABSlider.setUnit((Kelvin.UNIT));
		ABSlider.setShowBorder(true);
		ABSlider.setLabel("RB epsilon");
		ABSlider.setMinimum(eMin);
		ABSlider.setMaximum(eMax);
		ABSlider.setNMajor(5);
		ABSlider.getSlider().setSnapToTicks(true);
		ABSlider.setPostAction(fullResetAction);
		//    final DeviceSlider BBSlider = new DeviceSlider(BBbonded, "epsilon");
		final DeviceSlider BBSlider = new DeviceSlider(sim.getController(), new ModifierSQWEpsilon(sim.BBbonded));
		BBSlider.setUnit((Kelvin.UNIT));
		BBSlider.setShowBorder(true);
		BBSlider.setLabel("BB epsilon");
		BBSlider.setMinimum(eMin);
		BBSlider.setMaximum(eMax);
		BBSlider.setNMajor(5);
		BBSlider.getSlider().setSnapToTicks(true);
		BBSlider.setPostAction(fullResetAction);

		DiameterModifier sizeModifier = new DiameterModifier(sim.AAbonded,
				sim.ABbonded, sim.BBbonded, sim.speciesA, sim.speciesB);
		DeviceSlider sizeSlider = new DeviceSlider(sim.getController(), sizeModifier);
		//       sizeSlider.setShowBorder(true); //doesn't look right with just one
		// slider in
		// pane
		sizeSlider.setLabel("Atom size");
		sizeSlider.setPrecision(2);
		sizeSlider.setMinimum(0.0);
		sizeSlider.setMaximum(2.0);
		sizeSlider.setNMajor(4);
		sizeSlider.setValue(2.0);
		sizeSlider.setShowValues(true);
		sizeSlider.setEditValues(true);
		sizeSlider.setPostAction(fullResetAction);

		//sliders to adjust widths of wells
		eMin = 10;
		eMax = 100;
		int majorSpacing = 15;
		int minorSpacing = 5;
		DeviceSlider AAWellSlider = new DeviceSlider(sim.getController(), new WellModifier(sim.AAbonded));
		DeviceSlider ABWellSlider = new DeviceSlider(sim.getController(), new WellModifier(sim.ABbonded));
		DeviceSlider BBWellSlider = new DeviceSlider(sim.getController(), new WellModifier(sim.BBbonded));
		AAWellSlider.setShowBorder(true);
		ABWellSlider.setShowBorder(true);
		BBWellSlider.setShowBorder(true);
		AAWellSlider.setLabel("RR core");
		ABWellSlider.setLabel("RB core");
		BBWellSlider.setLabel("BB core");
		AAWellSlider.setMinimum(eMin);
		AAWellSlider.setMaximum(eMax);
		ABWellSlider.setMinimum(eMin);
		ABWellSlider.setMaximum(eMax);
		BBWellSlider.setMinimum(eMin);
		BBWellSlider.setMaximum(eMax);
		AAWellSlider.getSlider().setSnapToTicks(true);
		ABWellSlider.getSlider().setSnapToTicks(true);
		BBWellSlider.getSlider().setSnapToTicks(true);
		AAWellSlider.getSlider().setValue(50);
		ABWellSlider.getSlider().setValue(50);
		BBWellSlider.getSlider().setValue(50);
		AAWellSlider.getSlider().setMajorTickSpacing(majorSpacing);
		AAWellSlider.getSlider().setMinorTickSpacing(minorSpacing);
		ABWellSlider.getSlider().setMajorTickSpacing(majorSpacing);
		ABWellSlider.getSlider().setMinorTickSpacing(minorSpacing);
		BBWellSlider.getSlider().setMajorTickSpacing(majorSpacing);
		BBWellSlider.getSlider().setMinorTickSpacing(minorSpacing);
		AAWellSlider.getSlider().setLabelTable(
				AAWellSlider.getSlider().createStandardLabels(majorSpacing));
		ABWellSlider.getSlider().setLabelTable(
				ABWellSlider.getSlider().createStandardLabels(majorSpacing));
		BBWellSlider.getSlider().setLabelTable(
				BBWellSlider.getSlider().createStandardLabels(majorSpacing));
		AAWellSlider.setPostAction(fullResetAction);
		ABWellSlider.setPostAction(fullResetAction);
		BBWellSlider.setPostAction(fullResetAction);

		//so that display is updated when slider changes atom sizes
		sizeModifier.setDisplay(getDisplayBox(sim.box));

//		DisplayTextBox tBox = new DisplayTextBox(sim.thermometer.getDataInfo());
//		DataPump tPump = new DataPump (sim.thermometer, tBox);
//        sim.integratorHard1.addIntervalAction(tPump);
//        sim.integratorHard1.setActionInterval(tPump, 100);
//        tBox.setUnit(Kelvin.UNIT);
//		tBox.setLabel("Measured Temperature");
//		tBox.setLabelPosition(CompassDirection.NORTH);

		final AccumulatorAverageCollapsing tempAccum = new AccumulatorAverageCollapsing();
		final DataPump tPump = new DataPump(sim.thermometer, tempAccum);
		IntegratorListenerAction tPumpListener = new IntegratorListenerAction(tPump);
		sim.integratorHard1.getEventManager().addListener(tPumpListener);
		tPumpListener.setInterval(10);
		tempAccum.setPushInterval(10);
		tPump.setDataSink(tempAccum);
		final DisplayTextBoxesCAE tBox = new DisplayTextBoxesCAE();
		tempAccum.addDataSink(tBox,
				new AccumulatorAverage.StatType[]{tempAccum.MOST_RECENT,
						tempAccum.AVERAGE, tempAccum.ERROR});
		tBox.setLabel("Measured Temperature (K)");
		tBox.setUnit(Kelvin.UNIT);
		tBox.setLabelPosition(CompassDirection.NORTH);
		tBox.setLabelType(LabelType.BORDER);
		getController().getDataStreamPumps().add(tPump);

// NOTE : THE FOLLOWING IMPLEMENTATION IS CAUSING THE GRAPHIC
// TO BE UPDATED.  NEED TO REMOVE THE UPDATE OF THE GRAPHIC
// FROM THIS IMPLEMENTATION OF AVERAGE DISPLAY.
		//display of averages
		DataFork dimerFork = new DataFork();
		dimerPump = new DataPump(sim.meterDimerFraction, dimerFork);
		IntegratorListenerAction dimerPumpListener = new IntegratorListenerAction(dimerPump);
		sim.integratorHard1.getEventManager().addListener(dimerPumpListener);
		dimerPumpListener.setInterval(100);
		getController().getDataStreamPumps().add(dimerPump);

		DataGroupFilter filter1 = new DataGroupFilter(0);
		dimerFork.addDataSink(filter1);

		AccumulatorAverageFixed dimerFractionAccum = new AccumulatorAverageFixed(100);
		dimerFractionAccum.setPushInterval(10);
		filter1.setDataSink(dimerFractionAccum);

		DisplayTable table = new DisplayTable();
		dimerFractionAccum.addDataSink(table.getDataTable().makeDataSink(dimerFractionAccum.getDataInfo()),
				new AccumulatorAverage.StatType[]{dimerFractionAccum.AVERAGE, dimerFractionAccum.ERROR});


		table.setColumnHeader(new DataTag[]{dimerFractionAccum.getTag(dimerFractionAccum.AVERAGE)}, "Average");
		table.setColumnHeader(new DataTag[]{dimerFractionAccum.getTag(dimerFractionAccum.ERROR)}, "Error");
		table.setLabel("Fractions");

		DataSplitter splitter = new DataSplitter();
		DataGroupFilter filter2 = new DataGroupFilter(0);
		dimerFork.addDataSink(filter2);
		filter2.setDataSink(splitter);

		//display for history of mole fractions
		DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integratorHard1);
		DisplayPlotXChart plot = new DisplayPlotXChart();
		plot.setLabel("Composition");
		plot.setDoLegend(true);
//        int nData = sim.meterDimerFraction.getDataInfo().getLength();
//        DataTable.DataInfoTable dimerInfo = (DataTable.DataInfoTable)sim.meterDimerFraction.getDataInfo();
		int nData = filter1.getDataInfo().getLength();
		DataTable.DataInfoTable dimerInfo = (DataTable.DataInfoTable) filter1.getDataInfo();

		final AccumulatorHistory[] fractionHistory = new AccumulatorHistory[nData];
		for (int i = 0; i < nData; i++) {
			fractionHistory[i] = new AccumulatorHistory();
			fractionHistory[i].setTimeDataSource(timeCounter);

			splitter.setDataSink(i, fractionHistory[i]);
			fractionHistory[i].addDataSink(plot.getDataSet().makeDataSink());
			plot.setLegend(new DataTag[]{fractionHistory[i].getTag()}, dimerInfo.getRowHeader(i));
		}

		DataGroupFilter filter3 = new DataGroupFilter(1);
		dimerFork.addDataSink(filter3);
		densityAccum = new AccumulatorAverageCollapsing();
		densityAccum.setPushInterval(10);
		filter3.setDataSink(densityAccum);

		densityDisplay = new DisplayTextBoxesCAE();
		densityAccum.addDataSink(densityDisplay,
				new AccumulatorAverage.StatType[]{densityAccum.MOST_RECENT,
						densityAccum.AVERAGE,
						densityAccum.ERROR});
		densityDisplay.setLabel("Molecular density (" + Angstrom.UNIT.symbol() + "^-2)");
		dimerPump.actionPerformed();
		densityDisplay.putData(densityAccum.getData());
		densityDisplay.setLabelType(LabelType.BORDER);

//        filter3.setDataSink(new DataSinkConsole());
		DeviceDelaySlider delaySlider = new DeviceDelaySlider(sim.getController());

		//************* Lay out components ****************//

		//panel for the species editors
		JPanel speciesEditors = new JPanel(new GridLayout(0, 1));
		speciesEditors.add(AEditor);
		speciesEditors.add(BEditor);
		speciesEditors.setBorder(new TitledBorder(
				null, "Species Adjustment", TitledBorder.CENTER, TitledBorder.TOP));

		JPanel massEditors = new JPanel(new GridLayout(0, 1));
		massEditors.add(AMassEditor);
		massEditors.add(BMassEditor);

		//panel of well-depth sliders
		JPanel epsilonSliders = new JPanel(new GridLayout(0, 1));
		epsilonSliders.add(AASlider.graphic());
		epsilonSliders.add(ABSlider.graphic());
		epsilonSliders.add(BBSlider.graphic());

		//panel of well-width sliders
		JPanel lambdaSliders = new JPanel(new GridLayout(0, 1));
		lambdaSliders.add(AAWellSlider.graphic());
		lambdaSliders.add(ABWellSlider.graphic());
		lambdaSliders.add(BBWellSlider.graphic());

		//panel for size slider
		JPanel sizeSliders = new JPanel(new GridLayout(0, 1));
		sizeSliders.add(sizeSlider.graphic());

		//tabbed pane for both sets of sliders
		final JTabbedPane sliderPanel = new JTabbedPane();
		sliderPanel.setBorder(new TitledBorder(
				null, "Potential Adjustment", TitledBorder.CENTER, TitledBorder.TOP));
		sliderPanel.add("Well depth (K)", epsilonSliders);
		sliderPanel.add("Core size (%)", lambdaSliders);
		sliderPanel.add("Atom size (\u00C5)", sizeSliders);
		sliderPanel.addChangeListener(new javax.swing.event.ChangeListener() {
			public void stateChanged(javax.swing.event.ChangeEvent event) {
				sliderPanel.invalidate();
				sliderPanel.validate();
			}
		});
		sliderPanel.add("Species", speciesEditors);
		sliderPanel.add("Mass (Da)", massEditors);

		//top panel for control, temperature, potential adjustment
		add(temperatureSelect);
		getPanel().controlPanel.add(sliderPanel, vertGBC);
		getPanel().controlPanel.add(delaySlider.graphic(), vertGBC);
		add(plot);
		add(table);
		add(tBox);
		add(densityDisplay);

		getController().getReinitButton().setPostAction(new IAction() {
			public void actionPerformed() {
				for (int i = 0; i < fractionHistory.length; i++) {
					fractionHistory[i].reset();
				}

				tPump.actionPerformed();
				tBox.putData(tempAccum.getData());

				dimerPump.actionPerformed();
				densityDisplay.putData(densityAccum.getData());

				getDisplayBox(sim.box).graphic().repaint();
			}
		});


		table.setPrecision(6);

		initializing = false;
	}

	public static void main(String[] args) {

		ReactionEquilibrium sim = new ReactionEquilibrium();
		ReactionEquilibriumGraphic graphic = new ReactionEquilibriumGraphic(sim, sim.getSpace());
		SimulationGraphic.makeAndDisplayFrame(graphic.getPanel(), APP_NAME);
	}

	//=================================================================
	//panel containing species-editing devices

	public class MySpeciesEditor extends JPanel {

		//	public DeviceSlider nSlider;
		public DeviceNSelector nSlider;

		public MySpeciesEditor(SpeciesGeneral s, String label) {
			super();
			nSlider = new DeviceNSelector(sim.getController());
			nSlider.setResetAction(new SimulationRestart(sim));
			nSlider.setSpecies(s);
			nSlider.setBox(sim.box);
			//nSlider.setDisplayBox(DisplayBox1);
			nSlider.setMinimum(0);
			nSlider.setMaximum(60);
			nSlider.setPostAction(new IAction() {
				public void actionPerformed() {
					AtomLeafAgentManager agentManager = sim.getAgentManager();
					AtomIteratorLeafAtoms iter = new AtomIteratorLeafAtoms(sim.box);
					iter.reset();
					for (IAtom a = iter.nextAtom(); a != null; a = iter.nextAtom()) {
						//                      System.out.println(iter.peek().toString());
						agentManager.setAgent(a, null);
					}
					try {
						sim.integratorHard1.reset();
					} catch (ConfigurationOverlapException e) {
					}
					getDisplayBox(sim.box).repaint();

					resetAction.actionPerformed();
					//yay for a push data model
					dimerPump.actionPerformed();
					densityDisplay.putData(densityAccum.getData());
				}
			});
			setLayout(new FlowLayout());
			add(nSlider.graphic());
			setBorder(new TitledBorder(label));
		}
	} //end of MySpeciesEditor

	public static class MassEditor extends JPanel {

		public final DeviceBox mass;

		public MassEditor(Controller controller, SpeciesGeneral species, String label) {
			super();
			mass = new DeviceBox(controller);
			//listener for changes to mass textbox
			mass.setModifier(new ModifierGeneral(species.getLeafType().getElement(), "mass"));
			mass.setInteger(true);
			setLayout(new FlowLayout());
			add(mass.graphic());
			setBorder(new TitledBorder(label));
		}
	} //end of MySpeciesEditor

	//---------------------------------------------

	class DiameterModifier implements Modifier {
		P2SquareWellBonded potentialRR, potentialRB, potentialBB;

		SpeciesGeneral speciesR, speciesB;

		DisplayBox display;

		DiameterModifier(P2SquareWellBonded potentialRR,
                         P2SquareWellBonded potentialRB, P2SquareWellBonded potentialBB,
                         SpeciesGeneral speciesR, SpeciesGeneral speciesB) {
			this.potentialRR = potentialRR;
			this.potentialRB = potentialRB;
			this.potentialBB = potentialBB;
			this.speciesR = speciesR;
			this.speciesB = speciesB;
		}

		public Dimension getDimension() {
			return Length.DIMENSION;
		}

		public void setValue(double d) {
			if (d == 0.0)
				d = 0.01;
			double oldValue = potentialRR.getCollisionDiameter(1);
			double lambda = oldValue / potentialRR.getCollisionDiameter(0);
			potentialRR.setCollisionDiameter(0, d / lambda);
			lambda = oldValue / potentialRB.getCollisionDiameter(0);
			potentialRB.setCollisionDiameter(0, d / lambda);
			lambda = oldValue / potentialBB.getCollisionDiameter(0);
			potentialBB.setCollisionDiameter(0, d / lambda);
			potentialRR.setCollisionDiameter(1, d);
			potentialRB.setCollisionDiameter(1, d);
			potentialBB.setCollisionDiameter(1, d);
			sim.integratorHard1.setMaxCollisionDiameter(speciesB.getLeafType(), d);
			sim.integratorHard1.setMaxCollisionDiameter(speciesR.getLeafType(), d);
			((DiameterHashByType) getDisplayBox(sim.box).getDiameterHash()).setDiameter(speciesR.getLeafType(), d);
			((DiameterHashByType) getDisplayBox(sim.box).getDiameterHash()).setDiameter(speciesB.getLeafType(), d);
			sim.resetBonding();
			if (display != null)
				display.repaint();
		}

		public double getValue() {
			return ((DiameterHashByType) getDisplayBox(sim.box).getDiameterHash()).getDiameter(speciesR.getLeafType());
		}

		public void setDisplay(DisplayBox display) {
			this.display = display;
		}

		public String getLabel() {
			return "Diameter";
		}
	}

	class WellModifier implements Modifier {

		P2HardGeneric potential;

		WellModifier(P2HardGeneric pot) {
			potential = pot;
		}

		public String getLabel() {
			return "WellMod";
		}

		public Dimension getDimension() {
			return Null.DIMENSION;
		}

		public synchronized void setValue(double d) {
			if (initializing)
				return;
			potential.setCollisionDiameter(0, 0.01 * d * potential.getCollisionDiameter(1));
		}

		public double getValue() {
			return 100.0 * potential.getCollisionDiameter(0) / potential.getCollisionDiameter(1);
		}
	}//end of WellModulator

}
