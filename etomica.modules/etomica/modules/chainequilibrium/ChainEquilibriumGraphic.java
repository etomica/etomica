package etomica.modules.chainequilibrium;

import java.awt.GridBagConstraints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.SwingConstants;

import etomica.api.IAction;
import etomica.api.IAtomLeaf;
import etomica.api.IAtomSet;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.DataSinkTable;
import etomica.data.DataTag;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceDelaySlider;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTable;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.modifier.ModifierGeneral;
import etomica.space.ISpace;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.util.HistoryCollapsingAverage;
import etomica.util.Constants.CompassDirection;

/**
 * Module for chain reaction (polymerization) using ChainEquilibriumSim as the
 * simulation class.  Original module by William Scharmach and Matt Moynihan.
 * Later revamped based on module redesign by William M. Chirdon.
 * 
 * @author William Scharmach
 * @author Matt Moynihan
 * @author Andrew Schultz
 */
public class ChainEquilibriumGraphic extends SimulationGraphic {

	private static final String APP_NAME = "Chain Reaction Equilibrium";
	private static final int REPAINT_INTERVAL = 1;

    protected ChainEquilibriumSim sim;

    public ChainEquilibriumGraphic(ChainEquilibriumSim simulation, ISpace _space) {

		super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL, _space);
        this.sim = simulation;
        
        ArrayList dataStreamPumps = getController().getDataStreamPumps();
        
        getDisplayBox(sim.box).setPixelUnit(new Pixel(10));

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        // ********* Data Declaration Section *******	
        int eMin = 0, eMax = 4000;

        // **** Stuff that Modifies the Simulation

        final IAction resetAction = getController().getSimRestart().getDataResetAction();
        
        DeviceDelaySlider delaySlider = new DeviceDelaySlider(sim.controller1, sim.activityIntegrate);

        // Sliders on Well depth page
        final DeviceSlider ABSlider = sliders(eMin, eMax, "Diol-Carboxylic Acid", sim.ABbonded);
        final DeviceSlider ACSlider = sliders(eMin, eMax, "Diol-Crosslinker", sim.ACbonded);
        ABSlider.setPostAction(resetAction);
        ACSlider.setPostAction(resetAction);
        
        DisplayTextBox tBox = new DisplayTextBox();

        JPanel speciesEditors = new JPanel(new java.awt.GridLayout(0, 1));
        JPanel epsilonSliders = new JPanel(new java.awt.GridLayout(0, 1));

        DataPump tPump = new DataPump (sim.thermometer, tBox);
        add(tBox);
        dataStreamPumps.add(tPump);

        final MeterChainLength molecularCount = new MeterChainLength(sim.agentManager);
        molecularCount.setBox(sim.box);
        AccumulatorAverage accumulator = new AccumulatorAverageFixed(10);
        accumulator.setPushInterval(10);
        DataFork mwFork = new DataFork();
        DataPump pump = new DataPump(molecularCount,mwFork);
        mwFork.addDataSink(accumulator);
        dataStreamPumps.add(pump);
        sim.integratorHard.addIntervalAction(pump);
        sim.integratorHard.setActionInterval(pump, 10);
        
        MolecularWeightAvg molecularWeightAvg = new MolecularWeightAvg();
        mwFork.addDataSink(molecularWeightAvg);
        AccumulatorAverageCollapsing mwAvg = new AccumulatorAverageCollapsing();
        mwAvg.setPushInterval(1);
        molecularWeightAvg.setDataSink(mwAvg);

        MolecularWeightAvg2 molecularWeightAvg2 = new MolecularWeightAvg2();
        mwFork.addDataSink(molecularWeightAvg2);
        AccumulatorAverageCollapsing mwAvg2 = new AccumulatorAverageCollapsing();
        mwAvg2.setPushInterval(1);
        molecularWeightAvg2.setDataSink(mwAvg2);

        MonomerConversion monomerConversion = new MonomerConversion();
        mwFork.addDataSink(monomerConversion);
        AccumulatorHistory monomerConversionHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        monomerConversion.setDataSink(monomerConversionHistory);

        MeterConversion reactionConversion = new MeterConversion(sim.agentManager);
        AccumulatorHistory conversionHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        DataPump conversionPump = new DataPump(reactionConversion, conversionHistory);
        sim.integratorHard.addIntervalAction(conversionPump);
        sim.integratorHard.setActionInterval(conversionPump, 10);
        dataStreamPumps.add(conversionPump);

        sim.integratorHard.setTimeStep(0.01);

        DataSinkTable dataTable = new DataSinkTable();
        accumulator.addDataSink(dataTable.makeDataSink(),new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
        DisplayTable THING = new DisplayTable(dataTable);
        THING.setTransposed(false);
        THING.setShowingRowLabels(false);
        THING.setPrecision(7);

        DisplayPlot compositionPlot = new DisplayPlot();
        accumulator.addDataSink(compositionPlot.getDataSet().makeDataSink(),new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
        compositionPlot.setDoLegend(false);

        DisplayPlot conversionPlot = new DisplayPlot();
        monomerConversionHistory.addDataSink(conversionPlot.getDataSet().makeDataSink());
        conversionPlot.setLegend(new DataTag[]{monomerConversion.getTag()}, "monomer conversion");    
        conversionHistory.addDataSink(conversionPlot.getDataSet().makeDataSink());
        conversionPlot.setLegend(new DataTag[]{reactionConversion.getTag()}, "reaction conversion");    

        DisplayTextBoxesCAE mwBox = new DisplayTextBoxesCAE();
        mwBox.setAccumulator(mwAvg);
        mwBox.setLabel("Number Avg Molecular Weight");
        add(mwBox);

        DisplayTextBoxesCAE mw2Box = new DisplayTextBoxesCAE();
        mw2Box.setAccumulator(mwAvg2);
        mw2Box.setLabel("Weight Avg Molecular Weight");
        add(mw2Box);

        getController().getReinitButton().setPostAction(new IAction() {
            public void actionPerformed() {
                resetBonds();
                getDisplayBox(sim.box).repaint();
                molecularCount.reset();
            }
        });

        // Things to Do while simulation is on (It is adding the DataPumps, which run off the meters)
        sim.integratorHard.addIntervalAction(tPump);
        // Setting up how often it operates. 
        sim.integratorHard.setActionInterval(tPump, 10);

        DeviceThermoSlider temperatureSelect = new DeviceThermoSlider(sim.controller1);
        temperatureSelect.setUnit(Kelvin.UNIT);
        temperatureSelect.setIntegrator(sim.integratorHard);
        temperatureSelect.setTemperature(150);
        temperatureSelect.setMaximum(1200);
        temperatureSelect.setIsothermal();
        temperatureSelect.setSliderPostAction(resetAction);
        temperatureSelect.addRadioGroupActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                resetAction.actionPerformed();
            }
        });
        
        ColorSchemeByType colorScheme = (ColorSchemeByType)getDisplayBox(sim.box).getColorScheme();
        colorScheme.setColor(sim.speciesA.getLeafType(), java.awt.Color.RED);
        colorScheme.setColor(sim.speciesB.getLeafType(), java.awt.Color.BLACK);
        colorScheme.setColor(sim.speciesC.getLeafType(), java.awt.Color.GREEN);

        DeviceNSelector nSliderA = new DeviceNSelector(sim.getController());
        nSliderA.setSpecies(sim.speciesA);
        nSliderA.setBox(sim.box);
        nSliderA.setShowBorder(true);
        nSliderA.setLabel("Diol");
        nSliderA.setNMajor(4);
        IAction reset = new IAction() {
            public void actionPerformed() {
                getController().getSimRestart().actionPerformed();
                resetBonds();
                getDisplayBox(sim.box).repaint();
                molecularCount.reset();
            }
        };
        nSliderA.setResetAction(reset);
        nSliderA.setMaximum(400);
        DeviceNSelector nSliderB = new DeviceNSelector(sim.getController());
        nSliderB.setSpecies(sim.speciesB);
        nSliderB.setBox(sim.box);
        nSliderB.setShowBorder(true);
        nSliderB.setLabel("Carboxylic Acid");
        nSliderB.setResetAction(reset);
        nSliderB.setNMajor(4);
        nSliderB.setMaximum(400);
        DeviceNSelector nSliderC = new DeviceNSelector(sim.getController());
        nSliderC.setSpecies(sim.speciesC);
        nSliderC.setBox(sim.box);
        nSliderC.setShowBorder(true);
        nSliderC.setLabel("Crosslinker");
        nSliderC.setResetAction(reset);
        nSliderC.setNMajor(4);
        nSliderC.setMaximum(20);

        tBox.setUnit(Kelvin.UNIT);
        tBox.setLabel("Measured Temperature");
        tBox.setLabelPosition(CompassDirection.NORTH);

        compositionPlot.setLabel("Composition");
        conversionPlot.setLabel("Conversion");

        getPanel().tabbedPane.add(compositionPlot.getLabel(), compositionPlot.graphic());
        getPanel().tabbedPane.add(conversionPlot.getLabel(), conversionPlot.graphic());
        getPanel().tabbedPane.add("Averages", THING.graphic());


        speciesEditors.add(nSliderA.graphic());
        speciesEditors.add(nSliderB.graphic());
        speciesEditors.add(nSliderC.graphic());

        epsilonSliders.add(ABSlider.graphic(null));
        epsilonSliders.add(ACSlider.graphic(null));

        final JTabbedPane sliderPanel = new JTabbedPane();
        //panel for all the controls
        getPanel().controlPanel.add(delaySlider.graphic(), vertGBC);
        getPanel().controlPanel.add(temperatureSelect.graphic(), vertGBC);
        getPanel().controlPanel.add(sliderPanel, vertGBC);
        sliderPanel.add(epsilonSliders, "Well depth (K)");
        sliderPanel.add(speciesEditors, "Number of Molecules");

        //set the number of significant figures displayed on the table.
        javax.swing.table.DefaultTableCellRenderer numberRenderer = new javax.swing.table.DefaultTableCellRenderer() {
            java.text.NumberFormat formatter;
            {
                formatter = java.text.NumberFormat.getInstance();
                formatter.setMaximumFractionDigits(6);
            }

            public void setValue(Object value) {
                setText((value == null) ? "" : formatter.format(value));
            }
        };

        numberRenderer.setHorizontalAlignment(SwingConstants.RIGHT);
    }
    
    public void resetBonds() {
        IAtomSet atoms = sim.box.getLeafList();
        for (int i=0; i<atoms.getAtomCount(); i++) {
            IAtomLeaf[] bonds = (IAtomLeaf[])sim.agentManager.getAgent(atoms.getAtom(i));
            for (int j=0; j<bonds.length; j++) {
                bonds[j] = null;
            }
        }
    }

    public DeviceSlider sliders(int eMin, int eMax, String s, P2SquareWellBonded p){

        DeviceSlider AASlider = new DeviceSlider(sim.getController(), new ModifierGeneral(p, "epsilon"));
        AASlider.setUnit((Kelvin.UNIT));
        AASlider.setShowBorder(true);
        AASlider.setLabel(s);
        AASlider.setMinimum(eMin);
        AASlider.setMaximum(eMax);
        AASlider.setNMajor(4);
//        AASlider.getSlider().setSnapToTicks(true);

        return AASlider;
    }

    public static void main(String[] args) {
        ChainEquilibriumSim sim = new ChainEquilibriumSim();
        ChainEquilibriumGraphic graphic = new ChainEquilibriumGraphic(sim, sim.getSpace());
        SimulationGraphic.makeAndDisplayFrame(graphic.getPanel(), APP_NAME);
    }

    public static class Applet extends javax.swing.JApplet {

        public void init() {
			getRootPane().putClientProperty("defeatSystemEventQueueCheck", Boolean.TRUE);
	        ChainEquilibriumSim sim = new ChainEquilibriumSim();
			getContentPane().add(new ChainEquilibriumGraphic(sim, sim.getSpace()).getPanel());
        }
    }
}