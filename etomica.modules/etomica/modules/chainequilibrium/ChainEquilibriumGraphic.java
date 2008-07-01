package etomica.modules.chainequilibrium;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import javax.swing.BoxLayout;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.SwingConstants;

import etomica.action.SimulationRestart;
import etomica.api.IAction;
import etomica.api.IAtomLeaf;
import etomica.api.IAtomSet;
import etomica.api.ISpecies;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.DataSinkTable;
import etomica.data.DataSourceCountTime;
import etomica.data.DataTag;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceBox;
import etomica.graphics.DeviceDelaySlider;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTable;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.DisplayTimer;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.space.ISpace;
import etomica.units.Dimension;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.units.Mole;
import etomica.units.Pixel;
import etomica.units.Prefix;
import etomica.units.PrefixedUnit;
import etomica.units.Quantity;
import etomica.units.UnitRatio;
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
	private static final int REPAINT_INTERVAL = 2;

    protected ChainEquilibriumSim sim;

    public ChainEquilibriumGraphic(ChainEquilibriumSim simulation, ISpace _space) {

		super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL, _space);
        this.sim = simulation;
        
        int dataInterval = (int) (.04 / sim.integratorHard.getTimeStep());
        
        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();
        
        getDisplayBox(sim.box).setPixelUnit(new Pixel(10));

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        // ********* Data Declaration Section *******	
        int eMin = 0, eMax = 40;

        // **** Stuff that Modifies the Simulation

        final IAction resetAction = getController().getSimRestart().getDataResetAction();
        
        DeviceDelaySlider delaySlider = new DeviceDelaySlider(sim.controller1, sim.activityIntegrate);

        // Sliders on Well depth page
        final DeviceSlider ABSlider = sliders(eMin, eMax, "Diol-Carboxylic Acid", sim.ABbonded);
        final DeviceSlider ACSlider = sliders(eMin, eMax, "Diol-Crosslinker", sim.ACbonded);
        ABSlider.setPostAction(resetAction);
        ACSlider.setPostAction(resetAction);
        
        DeviceBox solventThermoFrac = new DeviceBox();
        solventThermoFrac.setController(sim.getController());
        solventThermoFrac.setModifier(new ModifierGeneral(new P2SquareWellBonded[]{sim.ABbonded, sim.ACbonded}, "solventThermoFrac"));
        solventThermoFrac.setLabel("fraction heat transfer to solvent");
        DisplayTextBox tBox = new DisplayTextBox();

        DisplayTimer displayTimer = new DisplayTimer(sim.integratorHard);
        add(displayTimer);
        
        DataSourceCountTime timer = new DataSourceCountTime(sim.integratorHard);

        DataFork tFork = new DataFork();
        final DataPump tPump = new DataPump (sim.thermometer, tFork);
        tFork.addDataSink(tBox);
        add(tBox);
        dataStreamPumps.add(tPump);
        AccumulatorHistory tHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        tHistory.setTimeDataSource(timer);
        tFork.addDataSink(tHistory);
        DisplayPlot tPlot = new DisplayPlot();
        tHistory.addDataSink(tPlot.getDataSet().makeDataSink());
        tPlot.setUnit(Kelvin.UNIT);
        tPlot.setLabel("Temperature");
        tPlot.getPlot().setYLabel("Temperature (K)");
        tPlot.setDoLegend(false);
        add(tPlot);
        sim.integratorHard.addIntervalAction(tPump);
        sim.integratorHard.setActionInterval(tPump, dataInterval);

        final MeterChainLength molecularCount = new MeterChainLength(sim.agentManager);
        molecularCount.setBox(sim.box);
        AccumulatorAverage accumulator = new AccumulatorAverageFixed(10);
        accumulator.setPushInterval(10);
        DataFork mwFork = new DataFork();
        final DataPump mwPump = new DataPump(molecularCount,mwFork);
        mwFork.addDataSink(accumulator);
        dataStreamPumps.add(mwPump);
        sim.integratorHard.addIntervalAction(mwPump);
        sim.integratorHard.setActionInterval(mwPump, dataInterval);
        
        MolecularWeightAvg molecularWeightAvg = new MolecularWeightAvg();
        mwFork.addDataSink(molecularWeightAvg);
        DataFork mwAvgFork = new DataFork();
        molecularWeightAvg.setDataSink(mwAvgFork);
        AccumulatorAverageCollapsing mwAvg = new AccumulatorAverageCollapsing();
        mwAvg.setPushInterval(1);
        mwAvgFork.addDataSink(mwAvg);
        final AccumulatorHistory mwHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        mwHistory.setTimeDataSource(timer);
        mwAvgFork.addDataSink(mwHistory);

        MolecularWeightAvg2 molecularWeightAvg2 = new MolecularWeightAvg2();
        mwFork.addDataSink(molecularWeightAvg2);
        DataFork mwAvg2Fork = new DataFork();
        molecularWeightAvg2.setDataSink(mwAvg2Fork);
        AccumulatorAverageCollapsing mwAvg2 = new AccumulatorAverageCollapsing();
        mwAvg2.setPushInterval(1);
        mwAvg2Fork.addDataSink(mwAvg2);
        final AccumulatorHistory mw2History = new AccumulatorHistory(new HistoryCollapsingAverage());
        mw2History.setTimeDataSource(timer);
        mwAvg2Fork.addDataSink(mw2History);

        MeterConversion reactionConversionDiol = new MeterConversion(sim.box, sim.agentManager);
        reactionConversionDiol.setSpecies(new ISpecies[]{sim.speciesA});
        final HistoryCollapsingAverage conversionHistoryDiol = new HistoryCollapsingAverage();
        AccumulatorHistory conversionHistoryAccDiol = new AccumulatorHistory(conversionHistoryDiol);
        conversionHistoryAccDiol.setTimeDataSource(timer);
        final DataPump conversionPumpDiol = new DataPump(reactionConversionDiol, conversionHistoryAccDiol);
        sim.integratorHard.addIntervalAction(conversionPumpDiol);
        sim.integratorHard.setActionInterval(conversionPumpDiol, dataInterval);
        dataStreamPumps.add(conversionPumpDiol);

        MeterConversion reactionConversionAcid = new MeterConversion(sim.box, sim.agentManager);
        reactionConversionAcid.setSpecies(new ISpecies[]{sim.speciesB, sim.speciesC});
        final HistoryCollapsingAverage conversionHistoryAcid = new HistoryCollapsingAverage();
        AccumulatorHistory conversionHistoryAccAcid = new AccumulatorHistory(conversionHistoryAcid);
        conversionHistoryAccAcid.setTimeDataSource(timer);
        final DataPump conversionPumpAcid = new DataPump(reactionConversionAcid, conversionHistoryAccAcid);
        sim.integratorHard.addIntervalAction(conversionPumpAcid);
        sim.integratorHard.setActionInterval(conversionPumpAcid, dataInterval);
        dataStreamPumps.add(conversionPumpAcid);

        final IAction resetData = new IAction() {
            public void actionPerformed() {
                sim.integratorHard.resetTime();
                molecularCount.reset();
                conversionPumpDiol.actionPerformed();
                conversionPumpAcid.actionPerformed();
                mwPump.actionPerformed();
                tPump.actionPerformed();
            }
        };

        getController().getResetAveragesButton().setLabel("Reset");
        getController().getResetAveragesButton().setPostAction(resetData);

        DataSinkTable dataTable = new DataSinkTable();
        accumulator.addDataSink(dataTable.makeDataSink(),new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
        DisplayTable THING = new DisplayTable(dataTable);
        THING.setTransposed(false);
        THING.setShowingRowLabels(false);
        THING.setPrecision(7);

        DisplayPlot compositionPlot = new DisplayPlot();
        accumulator.addDataSink(compositionPlot.getDataSet().makeDataSink(),new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
        compositionPlot.setDoLegend(false);

        DisplayPlot mwPlot = new DisplayPlot();
        mwPlot.setLabel("Molecular Weight");
        mwHistory.addDataSink(mwPlot.getDataSet().makeDataSink());
        mwPlot.setLegend(new DataTag[]{mwHistory.getTag()}, "Number Avg");
        mw2History.addDataSink(mwPlot.getDataSet().makeDataSink());
        mwPlot.setLegend(new DataTag[]{mw2History.getTag()}, "Weight Avg");

        DisplayPlot conversionPlot = new DisplayPlot();
        conversionHistoryAccDiol.addDataSink(conversionPlot.getDataSet().makeDataSink());
        conversionPlot.setLegend(new DataTag[]{reactionConversionDiol.getTag()}, "diol conversion");
        conversionHistoryAccAcid.addDataSink(conversionPlot.getDataSet().makeDataSink());
        conversionPlot.setLegend(new DataTag[]{reactionConversionAcid.getTag()}, "acid conversion");

        DisplayTextBoxesCAE mwBox = new DisplayTextBoxesCAE();
        mwBox.setAccumulator(mwAvg);
        mwBox.setLabel("Number Avg Molecular Weight");
        add(mwBox);

        DisplayTextBoxesCAE mw2Box = new DisplayTextBoxesCAE();
        mw2Box.setAccumulator(mwAvg2);
        mw2Box.setLabel("Weight Avg Molecular Weight");
        add(mw2Box);

        ((SimulationRestart)getController().getReinitButton().getAction()).setConfiguration(sim.config);
        getController().getReinitButton().setPostAction(new IAction() {
            public void actionPerformed() {
                resetBonds();
                getDisplayBox(sim.box).repaint();
                resetData.actionPerformed();
            }
        });

        DeviceThermoSlider temperatureSelect = new DeviceThermoSlider(sim.controller1);
        temperatureSelect.setUnit(Kelvin.UNIT);
        temperatureSelect.setIntegrator(sim.integratorHard);
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
        nSliderA.setLabel("Di-ols (red)");
        nSliderA.setNMajor(4);
        IAction reset = new IAction() {
            public void actionPerformed() {
                getController().getSimRestart().actionPerformed();
                resetBonds();
                getDisplayBox(sim.box).repaint();
                resetData.actionPerformed();
            }
        };
        nSliderA.setResetAction(reset);
        nSliderA.setMaximum(1000);
        nSliderA.setShowValues(true);
        nSliderA.setShowSlider(false);
        nSliderA.setEditValues(true);
        DeviceNSelector nSliderB = new DeviceNSelector(sim.getController());
        nSliderB.setSpecies(sim.speciesB);
        nSliderB.setBox(sim.box);
        nSliderB.setShowBorder(true);
        nSliderB.setLabel("Di-acid (black)");
        nSliderB.setResetAction(reset);
        nSliderB.setNMajor(4);
        nSliderB.setMaximum(1000);
        nSliderB.setShowValues(true);
        nSliderB.setShowSlider(false);
        nSliderB.setEditValues(true);
        DeviceNSelector nSliderC = new DeviceNSelector(sim.getController());
        nSliderC.setSpecies(sim.speciesC);
        nSliderC.setBox(sim.box);
        nSliderC.setShowBorder(true);
        nSliderC.setLabel("Crosslinker (green)");
        nSliderC.setResetAction(reset);
        nSliderC.setNMajor(4);
        nSliderC.setMaximum(50);
        nSliderC.setShowValues(true);
        nSliderC.setShowSlider(false);
        nSliderC.setEditValues(true);

        tBox.setUnit(Kelvin.UNIT);
        tBox.setLabel("Measured Temperature");
        tBox.setLabelPosition(CompassDirection.NORTH);

        compositionPlot.setLabel("Composition");
        conversionPlot.setLabel("Conversion");

        add(compositionPlot);
        add(mwPlot);
        JPanel conversionPanel = new JPanel(new GridBagLayout());
        conversionPanel.add(conversionPlot.graphic(), vertGBC);
        
        DeviceBox conversionHistoryLength = new DeviceBox();
        conversionHistoryLength.setInteger(true);
        conversionHistoryLength.setController(sim.getController());
        conversionHistoryLength.setModifier(new Modifier() {

            public Dimension getDimension() {
                return Quantity.DIMENSION;
            }

            public String getLabel() {
                return "history length";
            }

            public double getValue() {
                return conversionHistoryDiol.getHistoryLength();
            }

            public void setValue(double newValue) {
                conversionHistoryDiol.setHistoryLength((int)newValue);
                conversionHistoryAcid.setHistoryLength((int)newValue);
            }
        });
        conversionPanel.add(conversionHistoryLength.graphic(),vertGBC);
        
        getPanel().tabbedPane.add("Conversion" , conversionPanel);

        JPanel speciesEditors = new JPanel(new java.awt.GridLayout(0, 1));
        JPanel epsilonSliders = new JPanel(new java.awt.GridBagLayout());

        speciesEditors.add(nSliderA.graphic());
        speciesEditors.add(nSliderB.graphic());
        speciesEditors.add(nSliderC.graphic());

        epsilonSliders.add(ABSlider.graphic(null), vertGBC);
        epsilonSliders.add(ACSlider.graphic(null), vertGBC);
        epsilonSliders.add(solventThermoFrac.graphic(), vertGBC);

        final JTabbedPane sliderPanel = new JTabbedPane();
        //panel for all the controls
        getPanel().controlPanel.add(delaySlider.graphic(), vertGBC);
        getPanel().controlPanel.add(temperatureSelect.graphic(), vertGBC);
        getPanel().controlPanel.add(sliderPanel, vertGBC);
        sliderPanel.add(epsilonSliders, "Reaction Energy (kJ/mol)");
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
        AASlider.setUnit(new UnitRatio(new PrefixedUnit(Prefix.KILO, Joule.UNIT), Mole.UNIT));
        AASlider.doUpdate();
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