/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.statistics;

import etomica.action.IAction;
import etomica.action.ResetAccumulators;
import etomica.action.SimulationDataAction;
import etomica.data.*;
import etomica.data.histogram.HistogramCollapsing;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.meter.*;
import etomica.graphics.*;
import etomica.modifier.ModifierBoolean;
import etomica.modules.ensembles.LJMC;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;

public class StatisticsMCGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Statistics";
    private final static int REPAINT_INTERVAL = 100;
    protected final LJMC sim;

    protected boolean volumeChanges = false;
    protected boolean constMu = false;

    public StatisticsMCGraphic(final LJMC simulation, int moduleNum) {

        super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL, simulation.getSpace(), simulation.getController());

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

        this.sim = simulation;

        //display of box, timer
        ColorSchemeByType colorScheme = new ColorSchemeByType(sim);
        colorScheme.setColor(sim.species.getLeafType(), Color.red);
        getDisplayBox(sim.box).setColorScheme(new ColorSchemeByType(sim));

        DataSourceCountSteps timeCounter = new DataSourceCountSteps(sim.integrator);

        // Number density box
        JPanel historyPanel = null;
        java.awt.Dimension d = new Dimension(600, 400);

        final HistoryPlotBits dHPB, peHPB, pHPB, widomHPB;
        final DataPumpListener pPump, dPump, pePump;
        if (moduleNum == 1) {
            MeterDensity densityMeter = new MeterDensity(sim.getSpace());
            densityMeter.setBox(sim.box);
            dPump = new DataPumpListener(densityMeter, null, 100);
            dataStreamPumps.add(dPump);

            historyPanel = new JPanel(new GridLayout(0, 1));
            dHPB = makeHistoryPlot(dataStreamPumps, timeCounter, historyPanel, dPump, "Density");
            DisplayPlot dPlot = dHPB.plot;
            dHPB.avg.setPushInterval(10);

            MeterPotentialEnergyFromIntegrator peMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
            AccumulatorHistory peHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
            peHistory.setTimeDataSource(timeCounter);
            pePump = new DataPumpListener(peMeter, null, 1);

            peHPB = makeHistoryPlot(dataStreamPumps, timeCounter, historyPanel, pePump, "Potential Energy");
            DataFork peFork = peHPB.fork;
            peHPB.avg.setPushInterval(10);

            JScrollPane historyPane = new JScrollPane(historyPanel);

            // Add plots page to tabbed pane
            addAsTab(historyPane, "History", true);

            // Set the size of the plots and the scoll pane containing the plots.
            // Want 2 of the 3 plots displayed
            d = dPlot.getPlot().getPreferredSize();

            d.width += 40;
            d.height = d.height * 2 + 40;
            historyPane.setPreferredSize(d);

            MeterPressure pMeter = new MeterPressure(space);
            pMeter.setIntegrator(sim.integrator);
            pMeter.setBox(sim.box);
            pPump = new DataPumpListener(pMeter, null, 1000);
            pHPB = makeHistoryPlot(dataStreamPumps, timeCounter, historyPanel, pPump, "Pressure");

            addAsTab(createStatPanel(peFork, d, null, true), "Potential Energy", true);
        } else {
            dPump = pPump = pePump = null;
            dHPB = pHPB = peHPB = null;
        }

        MeterWidomInsertion meterWidom = new MeterWidomInsertion(space, sim.getRandom());
        meterWidom.setNInsert(1);
        meterWidom.setSpecies(sim.species);
        meterWidom.setResidual(false);
        meterWidom.setEnergyMeter(new MeterPotentialEnergy(sim.integrator.getPotentialMaster()));
        meterWidom.setBox(sim.box);
        meterWidom.setTemperature(sim.integrator.getTemperature());
        final DataPumpListener widomPump = new DataPumpListener(meterWidom, null, 1);
        DataFork widomFork = null;
        if (moduleNum == 1) {
            widomHPB = makeMuHistoryPlot(dataStreamPumps, timeCounter, historyPanel, widomPump, "Chemical Potential");
            widomFork = widomHPB.fork;
        } else {
            widomHPB = null;
            widomFork = new DataFork();
            widomPump.setDataSink(widomFork);
            sim.integrator.getEventManager().addListener(widomPump);
            dataStreamPumps.add(widomPump);
        }
        dataStreamPumps.add(widomPump);
        AccumulatorFactory muFactory = new AccumulatorFactory() {
            @Override
            public AccumulatorAverageFixed makeAccumulator() {
                return new AccumulatorMimicMu(sim.integrator);
            }

            public DataProcessor makeDataProcessor() {
                return new DataProcessorMu(sim.integrator);
            }
        };
        addAsTab(createHistogramPanel(widomFork, d, muFactory), "Chemical Potential", true);
//        AccumulatorAverageCollapsing widomAvg = new AccumulatorAverageCollapsing();
//        widomHPB.fork.addDataSink(widomAvg);
//        AccumulatorMimicMu accMu = new AccumulatorMimicMu(sim.integrator);
//        widomAvg.addDataSink(accMu);

        final IAction resetDisplays;
        if (moduleNum == 1) {
            final DisplayTextBoxesCAE dDisplay = new DisplayTextBoxesCAE();
            dDisplay.setAccumulator(dHPB.avg);
            final DisplayTextBoxesCAE pDisplay = new DisplayTextBoxesCAE();
            pDisplay.setAccumulator(pHPB.avg);
            final DisplayTextBoxesCAE peDisplay = new DisplayTextBoxesCAE();
            peDisplay.setAccumulator(peHPB.avg);
            final DisplayTextBoxesCAE muDisplay = new DisplayTextBoxesCAE();
            muDisplay.setAccumulator(widomHPB.avg);
            muDisplay.setDoShowCurrent(false);
            resetDisplays = new IAction() {
                public void actionPerformed() {
                    dPump.actionPerformed();
                    ;
                    dDisplay.putData(dHPB.avg.getData());
                    dDisplay.repaint();

                    pPump.actionPerformed();
                    pDisplay.putData(pHPB.avg.getData());
                    pDisplay.repaint();
                    pePump.actionPerformed();
                    peDisplay.putData(peHPB.avg.getData());
                    peDisplay.repaint();
                }
            };
            add(dDisplay);
            add(pDisplay);
            add(peDisplay);
            add(muDisplay);
        } else {
            resetDisplays = null;
        }

        //************* Lay out components ****************//

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        //temperature selector
        DeviceThermoSlider temperatureSelect = new DeviceThermoSlider(sim.getController(), sim.integrator);
        temperatureSelect.setPrecision(1);
        temperatureSelect.setMinimum(0.0);
        temperatureSelect.setMaximum(10.0);
        temperatureSelect.setSliderMajorValues(4);
        temperatureSelect.setIsothermalButtonsVisibility(false);

        IAction resetAction = null;
        if (moduleNum == 1) {
            resetAction = new IAction() {
                public void actionPerformed() {
                    meterWidom.setTemperature(sim.integrator.getTemperature());

                    resetDisplays.actionPerformed();
                    getDisplayBox(sim.box).graphic().repaint();
                }
            };
        } else {
            resetAction = new IAction() {
                public void actionPerformed() {
                    meterWidom.setTemperature(sim.integrator.getTemperature());
                    getDisplayBox(sim.box).graphic().repaint();
                }
            };
        }

        temperatureSelect.setSliderPostAction(resetAction);
        this.getController().getReinitButton().setPostAction(resetAction);
        this.getController().getResetAveragesButton().setPostAction(resetAction);

        getPanel().controlPanel.add(temperatureSelect.graphic(), vertGBC);

        JPanel pPanel = new JPanel(new GridBagLayout());

        final DeviceSlider pSlider = new DeviceSlider(sim.getController(), sim.mcMoveVolume, "pressure");
        pSlider.setEnabled(volumeChanges);
        pSlider.setMaximum(10);
        pSlider.setPrecision(2);
        pSlider.setNMajor(4);
        pSlider.setShowValues(true);
        pSlider.setEditValues(true);
        pSlider.setShowBorder(true);
        pSlider.setLabel("Pressure");

        DeviceCheckBox pCheckbox = new DeviceCheckBox("Volume changes", new ModifierBoolean() {

            public void setBoolean(boolean b) {
                if (b == volumeChanges) return;
                if (b) {
                    sim.integrator.getMoveManager().addMCMove(sim.mcMoveVolume);
                    pSlider.setEnabled(true);
                } else {
                    sim.integrator.getMoveManager().removeMCMove(sim.mcMoveVolume);
                    pSlider.setEnabled(false);
                }
                volumeChanges = b;
            }

            public boolean getBoolean() {
                return volumeChanges;
            }
        });
        pCheckbox.setController(sim.getController());

        pPanel.add(pCheckbox.graphic(), vertGBC);
        pPanel.add(pSlider.graphic(), vertGBC);
        getPanel().controlPanel.add(pPanel, vertGBC);


        JPanel muPanel = new JPanel(new GridBagLayout());

        final DeviceSlider muSlider = new DeviceSlider(sim.getController(), sim.mcMoveID, "mu");
        muSlider.setEnabled(constMu);
        muSlider.setMinimum(-10);
        muSlider.setMaximum(10);
        muSlider.setPrecision(1);
        muSlider.setNMajor(4);
        muSlider.setShowValues(true);
        muSlider.setEditValues(true);
        muSlider.setShowBorder(true);
        muSlider.setLabel("Chemical Potential");

        DeviceCheckBox muCheckbox = new DeviceCheckBox("Insert/Delete", new ModifierBoolean() {

            public void setBoolean(boolean b) {
                if (b == constMu) return;
                if (b) {
                    sim.integrator.getMoveManager().addMCMove(sim.mcMoveID);
                    muSlider.setEnabled(true);
                } else {
                    sim.integrator.getMoveManager().removeMCMove(sim.mcMoveID);
                    muSlider.setEnabled(false);
                }
                constMu = b;
            }

            public boolean getBoolean() {
                return constMu;
            }
        });
        muCheckbox.setController(sim.getController());

        muPanel.add(muCheckbox.graphic(), vertGBC);
        muPanel.add(muSlider.graphic(), vertGBC);
        getPanel().controlPanel.add(muPanel, vertGBC);

        DisplayTextBox vBox = new DisplayTextBox();
        vBox.setLabel("Volume");
        MeterVolume meterVolume = new MeterVolume();
        meterVolume.setBox(sim.box);
        DataPumpListener vPump = new DataPumpListener(meterVolume, vBox, 100);
        sim.integrator.getEventManager().addListener(vPump);

        DisplayTextBox nBox = new DisplayTextBox();
        nBox.setLabel("N");
        MeterNMolecules meterN = new MeterNMolecules();
        meterN.setBox(sim.box);
        DataPumpListener nPump = new DataPumpListener(meterN, nBox, 100);
        sim.integrator.getEventManager().addListener(nPump);

        JPanel vnPanel = new JPanel(new GridBagLayout());
        GridBagConstraints horizGBC = SimulationPanel.getHorizGBC();
        vnPanel.add(vBox.graphic(), horizGBC);
        vnPanel.add(nBox.graphic(), horizGBC);
        getPanel().controlPanel.add(vnPanel, vertGBC);

        final DeviceButton slowButton = new DeviceButton(sim.getController(), null);
        slowButton.setAction(new IAction() {
            public void actionPerformed() {
                int sleep = sim.activityIntegrate.getSleepPeriod();
                sleep = 1 - sleep;
                sim.activityIntegrate.setSleepPeriod(sleep);
                slowButton.setLabel(sleep == 0 ? "Slow" : "Fast");
            }
        });
        slowButton.setLabel("Slow");
        add(slowButton);

        ((SimulationDataAction) getController().getResetAveragesButton().getAction()).setStreamAction(new ResetAccumulators());
    }

    public static class HistoryPlotBits {
        public DisplayPlot plot;
        public DataFork fork;
        public AccumulatorHistory history;
        public AccumulatorAverage avg;
    }

    protected HistoryPlotBits makeHistoryPlot(ArrayList<DataPump> dataStreamPumps, DataSourceCountSteps timeCounter, JPanel historyPanel, DataPumpListener pump, String name) {
        HistoryPlotBits rv = new HistoryPlotBits();
        rv.history = new AccumulatorHistory(new HistoryCollapsingAverage());
        rv.history.setTimeDataSource(timeCounter);
        rv.avg = new AccumulatorAverageCollapsing(100);
        rv.fork = new DataFork(new IDataSink[]{rv.history, rv.avg});
        pump.setDataSink(rv.fork);
        sim.integrator.getEventManager().addListener(pump);
        rv.avg.setPushInterval(1);
        dataStreamPumps.add(pump);

        rv.plot = new DisplayPlot();
        rv.history.setDataSink(rv.plot.getDataSet().makeDataSink());
        rv.plot.setDoLegend(false);
        rv.plot.getPlot().setYLabel(name);
        rv.plot.setLabel(name);
        historyPanel.add(rv.plot.graphic());
        AccumulatorHistory[] pSinks = new AccumulatorHistory[3];
        for (int i = 0; i < pSinks.length; i++) {
            pSinks[i] = new AccumulatorHistory(new HistoryCollapsingDiscard());
            pSinks[i].setTimeDataSource(timeCounter);
        }
        DataProcessorBounds dpBounds = new DataProcessorBounds(pSinks);
        rv.avg.addDataSink(dpBounds, new AccumulatorAverage.StatType[]{rv.avg.AVERAGE, rv.avg.ERROR});
        String[] labels = new String[]{"avg-", "avg", "avg+"};
        for (int i = 0; i < pSinks.length; i++) {
            pSinks[i].addDataSink(rv.plot.getDataSet().makeDataSink());
            rv.plot.setLegend(new DataTag[]{pSinks[i].getTag()}, labels[i]);
        }
        return rv;
    }

    protected HistoryPlotBits makeMuHistoryPlot(ArrayList<DataPump> dataStreamPumps, DataSourceCountSteps timeCounter, JPanel historyPanel, DataPumpListener pump, String name) {
        HistoryPlotBits rv = new HistoryPlotBits();
        AccumulatorAverageCollapsing avg1 = new AccumulatorAverageCollapsing(100);
        AccumulatorMimicMu accMu = new AccumulatorMimicMu(sim.integrator);
        avg1.addDataSink(accMu);
        rv.avg = accMu;
        rv.fork = new DataFork(new IDataSink[]{rv.history, avg1});
        pump.setDataSink(rv.fork);
        sim.integrator.getEventManager().addListener(pump);
        accMu.setPushInterval(1);
        rv.avg.setPushInterval(1);
        dataStreamPumps.add(pump);

        rv.plot = new DisplayPlot();
        rv.plot.setDoLegend(false);
        rv.plot.getPlot().setYLabel(name);
        rv.plot.setLabel(name);
        historyPanel.add(rv.plot.graphic());
        AccumulatorHistory[] pSinks = new AccumulatorHistory[3];
        for (int i = 0; i < pSinks.length; i++) {
            pSinks[i] = new AccumulatorHistory(new HistoryCollapsingDiscard());
            pSinks[i].setTimeDataSource(timeCounter);
        }
        DataProcessorBounds dpBounds = new DataProcessorBounds(pSinks);
        rv.avg.addDataSink(dpBounds, new AccumulatorAverage.StatType[]{rv.avg.AVERAGE, rv.avg.ERROR});
        String[] labels = new String[]{"avg-", "avg", "avg+"};
        for (int i = 0; i < pSinks.length; i++) {
            pSinks[i].addDataSink(rv.plot.getDataSet().makeDataSink());
            rv.plot.setLegend(new DataTag[]{pSinks[i].getTag()}, labels[i]);
        }
        return rv;
    }

    protected JScrollPane createStatPanel(DataFork fork, Dimension paneSize, AccumulatorFactory accFactory, boolean doHistory) {
        JPanel panel = new JPanel(new GridLayout(0, 1));
        JScrollPane pane = new JScrollPane(panel);
        pane.setPreferredSize(paneSize);

        DataCollector collectorErr = new DataCollector();
        DataCollector collectorCor = new DataCollector();
        DataCollector collectorDifficulty = new DataCollector();
        DataCollector collectorSamples = new DataCollector();
        DataCollector collectorErrCorrected = new DataCollector();
        DataCollector collectorDiffCorrected = new DataCollector();
        DisplayPlot blockHistoryPlot = doHistory ? new DisplayPlot() : null;
        if (doHistory) {
            blockHistoryPlot.getPlot().setYLabel("Block Averages");
            panel.add(blockHistoryPlot.graphic());
            blockHistoryPlot.getDataSet().setUpdatingOnAnyChange(true);
        }
        for (int i = 0; i < 30; i++) {
            AccumulatorAverageFixed acc = new AccumulatorAverageFixed(1L << i);
            fork.addDataSink(acc);
            AccumulatorAverageFixed bit = acc;
            if (accFactory != null) {
                bit = accFactory.makeAccumulator();
                acc.addDataSink(bit);
            }
            bit.addDataSink(new DataAccPusher(i, collectorErr), new AccumulatorAverage.StatType[]{acc.ERROR});
            bit.addDataSink(new DataAccPusher(i, collectorCor), new AccumulatorAverage.StatType[]{acc.BLOCK_CORRELATION});
            bit.addDataSink(new DataAccDifficultyPusher(i, collectorDifficulty, acc), new AccumulatorAverage.StatType[]{acc.ERROR});
            bit.addDataSink(new DataAccSamplesPusher(i, collectorSamples, acc), new AccumulatorAverage.StatType[]{acc.STANDARD_DEVIATION, acc.ERROR});
            bit.addDataSink(new DataAccCorrectedPusher(i, collectorErrCorrected), new AccumulatorAverage.StatType[]{acc.ERROR, acc.BLOCK_CORRELATION});
            bit.addDataSink(new DataAccDiffCorrectedPusher(i, collectorDiffCorrected, acc), new AccumulatorAverage.StatType[]{acc.ERROR, acc.BLOCK_CORRELATION});
            if (i % 5 == 0) {
                DataFork blockFork = new DataFork();
                if (accFactory != null && false) {
                    DataProcessor dp = accFactory.makeDataProcessor();
                    acc.setBlockDataSink(dp);
                    dp.setDataSink(blockFork);
                } else {
                    acc.setBlockDataSink(blockFork);
                }

                if (doHistory) {
                    AccumulatorHistory accBlockHistory = new AccumulatorHistory(new HistoryStatistics(100));
                    blockFork.addDataSink(accBlockHistory);
                    accBlockHistory.addDataSink(blockHistoryPlot.getDataSet().makeDataSink());
                    accBlockHistory.setPushInterval(1 + 2000 / (1L << i));
                    blockHistoryPlot.setLegend(new DataTag[]{accBlockHistory.getTag()}, "" + (1L << i));
                }
            }
        }
        DisplayPlot peErrorPlot = new DisplayPlot();
        DataPumpListener peAccPumpErr = new DataPumpListener(collectorErr, peErrorPlot.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(peAccPumpErr);
        peErrorPlot.setLabel("PE Error");
        peErrorPlot.setLegend(new DataTag[]{collectorErr.getTag()}, "Computed");
        peErrorPlot.getPlot().setXLog(true);
        peErrorPlot.getPlot().setYLog(true);
        peErrorPlot.getPlot().setYLabel("Error");
        panel.add(peErrorPlot.graphic());
        DisplayPlot peCorPlot = new DisplayPlot();
        DataPumpListener peAccPumpCor = new DataPumpListener(collectorCor, peCorPlot.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(peAccPumpCor);
        peCorPlot.setLabel("PE Correlation");
        peCorPlot.setDoLegend(false);
        peCorPlot.getPlot().setYLabel("Correlation");
        peCorPlot.getPlot().setXLog(true);
        peCorPlot.getPlot().setYRange(-1, 1);
        panel.add(peCorPlot.graphic());
        DisplayPlot peDifficultyPlot = new DisplayPlot();
        DataPumpListener peAccPumpDifficulty = new DataPumpListener(collectorDifficulty, peDifficultyPlot.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(peAccPumpDifficulty);
        peDifficultyPlot.setLabel("PE Difficulty");
        peDifficultyPlot.setLegend(new DataTag[]{collectorDifficulty.getTag()}, "Computed");
        peDifficultyPlot.getPlot().setYLabel("Difficulty");
        peDifficultyPlot.getPlot().setXLog(true);
        peDifficultyPlot.getPlot().setYLog(true);
        panel.add(peDifficultyPlot.graphic());
        DisplayPlot peSamplesPlot = new DisplayPlot();
        DataPumpListener peAccPumpSamples = new DataPumpListener(collectorSamples, peSamplesPlot.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(peAccPumpSamples);
        peSamplesPlot.setLabel("PE Samples");
        peSamplesPlot.setDoLegend(false);
        peSamplesPlot.getPlot().setYLabel("Steps per Sample");
        peSamplesPlot.getPlot().setXLog(true);
        peSamplesPlot.getPlot().setYLog(true);
        panel.add(peSamplesPlot.graphic());
        DataPumpListener peAccPumpErrCorrected = new DataPumpListener(collectorErrCorrected, peErrorPlot.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(peAccPumpErrCorrected);
        peErrorPlot.setLegend(new DataTag[]{collectorErrCorrected.getTag()}, "Corrected");
        DataPumpListener peAccPumpDiffCorrected = new DataPumpListener(collectorDiffCorrected, peDifficultyPlot.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(peAccPumpDiffCorrected);
        peDifficultyPlot.setLegend(new DataTag[]{collectorDiffCorrected.getTag()}, "Corrected");

        return pane;
    }

    protected JScrollPane createHistogramPanel(DataFork fork, Dimension paneSize, AccumulatorFactory accFactory) {
        JPanel panel = new JPanel(new GridLayout(0, 1));
        JScrollPane pane = new JScrollPane(panel);
        pane.setPreferredSize(paneSize);

        DataCollector collectorErr = new DataCollector();
        DataCollector collectorCor = new DataCollector();
        DataCollector collectorDifficulty = new DataCollector();
        DataCollector collectorSamples = new DataCollector();
        DataCollector collectorErrCorrected = new DataCollector();
        DataCollector collectorDiffCorrected = new DataCollector();
        DisplayPlot blockHistogramPlot = new DisplayPlot();
        blockHistogramPlot.getPlot().setYLabel("Block Histogram");
        panel.add(blockHistogramPlot.graphic());
        blockHistogramPlot.getDataSet().setUpdatingOnAnyChange(true);
        if (accFactory != null) {
            blockHistogramPlot.getPlot().setYLog(true);
            blockHistogramPlot.getPlot().setXLog(true);
        }
        for (int i = 0; i < 30; i += 5) {
            AccumulatorAverageFixed acc = new AccumulatorAverageFixed(1L << i);
            fork.addDataSink(acc);
            AccumulatorAverageFixed bit = acc;
            if (accFactory != null) {
                bit = accFactory.makeAccumulator();
                acc.addDataSink(bit);
            }
            DataFork blockFork = new DataFork();
            if (accFactory != null && false) {
                DataProcessor dp = accFactory.makeDataProcessor();
                acc.setBlockDataSink(dp);
                dp.setDataSink(blockFork);
            } else {
                acc.setBlockDataSink(blockFork);
            }

            int nBins = 100 << ((30 - i) / 5);
            AccumulatorHistogram accBlockHistogram = new AccumulatorHistogram(new HistogramCollapsing(nBins), nBins);
            blockFork.addDataSink(accBlockHistogram);
            accBlockHistogram.addDataSink(blockHistogramPlot.getDataSet().makeDataSink());
            accBlockHistogram.setPushInterval(1 + 2000 / (1L << i));
            blockHistogramPlot.setLegend(new DataTag[]{accBlockHistogram.getTag()}, "" + (1L << i));
        }
        return pane;
    }

    public static void main(String[] args) {
        StatsParams params = new StatsParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        int D = params.D;
        int moduleNum = params.moduleNum;
        Space sp;
        if (D == 2) {
            sp = Space2D.getInstance();
        } else {
            sp = Space3D.getInstance();
        }

        StatisticsMCGraphic ljmcGraphic = new StatisticsMCGraphic(new LJMC(sp), moduleNum);
        SimulationGraphic.makeAndDisplayFrame(ljmcGraphic.getPanel(), APP_NAME);
    }

    public static class StatsParams extends ParameterBase {
        public int D = 3;
        public int moduleNum = 1;
    }
}
