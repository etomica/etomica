/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.statistics;

import etomica.action.IAction;
import etomica.action.ResetAccumulators;
import etomica.action.SimulationDataAction;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.meter.*;
import etomica.graphics.*;
import etomica.modifier.ModifierBoolean;
import etomica.modules.ensembles.LJMC;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;

public class StatisticsMCGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Statistics";
    private final static int REPAINT_INTERVAL = 100;
    protected final LJMC sim;

    protected boolean volumeChanges = false;
    protected boolean constMu = false;

    public StatisticsMCGraphic(final LJMC simulation, Space _space) {

        super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL, _space, simulation.getController());

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

        this.sim = simulation;

        //display of box, timer
        ColorSchemeByType colorScheme = new ColorSchemeByType(sim);
        colorScheme.setColor(sim.species.getLeafType(), Color.red);
        getDisplayBox(sim.box).setColorScheme(new ColorSchemeByType(sim));

        DataSourceCountSteps timeCounter = new DataSourceCountSteps(sim.integrator);

        // Number density box
        final MeterDensity densityMeter = new MeterDensity(sim.getSpace());
        densityMeter.setBox(sim.box);
        AccumulatorHistory dHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        dHistory.setTimeDataSource(timeCounter);
        final AccumulatorAverageCollapsing dAccumulator = new AccumulatorAverageCollapsing();
        dAccumulator.setPushInterval(10);
        double dmin = ((int) (densityMeter.getDataAsScalar() * 100)) * 0.01;
        if (dmin > 0) dmin -= 0.005;
        DataFork dFork = new DataFork(new IDataSink[]{dHistory, dAccumulator});
        final DataPumpListener dPump = new DataPumpListener(densityMeter, dFork, 100);
        sim.integrator.getEventManager().addListener(dPump);
        dHistory.setPushInterval(1);
        dataStreamPumps.add(dPump);

        DisplayPlot dPlot = new DisplayPlot();
        dHistory.setDataSink(dPlot.getDataSet().makeDataSink());
        dPlot.setDoLegend(false);
        dPlot.getPlot().setYLabel("Density");
        dPlot.setLabel("Density");

        MeterPotentialEnergyFromIntegrator peMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorHistory peHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        peHistory.setTimeDataSource(timeCounter);
        final AccumulatorAverageCollapsing peAccumulator = new AccumulatorAverageCollapsing();
        peAccumulator.setPushInterval(10);
        final DataPumpListener pePump = new DataPumpListener(peMeter, null, 1);

        JPanel historyPanel = new JPanel(new GridLayout(0, 1));
        makeHistoryPlot(dataStreamPumps, timeCounter, historyPanel, peAccumulator, pePump, "Pressure");

        historyPanel.add(dPlot.graphic());
        JScrollPane historyPane = new JScrollPane(historyPanel);

        // Add plots page to tabbed pane
        addAsTab(historyPane, "History", true);

        // Set the size of the plots and the scoll pane containing the plots.
        // Want 2 of the 3 plots displayed
        java.awt.Dimension d = dPlot.getPlot().getPreferredSize();

        d.width += 40;
        d.height = d.height * 2 + 40;
        historyPane.setPreferredSize(d);

        MeterPressure pMeter = new MeterPressure(space);
        pMeter.setIntegrator(sim.integrator);
        pMeter.setBox(sim.box);
        final AccumulatorAverageCollapsing pAccumulator = new AccumulatorAverageCollapsing();
        final DataPumpListener pPump = new DataPumpListener(pMeter, null, 1000);
        DataFork peFork = makeHistoryPlot(dataStreamPumps, timeCounter, historyPanel, pAccumulator, pPump, "Potential Energy");

        addAsTab(createStatPanel(peFork, d, null, true), "Potential Energy", true);

        MeterWidomInsertion meterWidom = new MeterWidomInsertion(space, sim.getRandom());
        meterWidom.setNInsert(1);
        meterWidom.setSpecies(sim.species);
        meterWidom.setResidual(false);
        meterWidom.setEnergyMeter(new MeterPotentialEnergy(sim.integrator.getPotentialMaster()));
        meterWidom.setBox(sim.box);
        meterWidom.setTemperature(sim.integrator.getTemperature());
        DataFork widomFork = new DataFork();
        final DataPumpListener widomPump = new DataPumpListener(meterWidom, widomFork, 1);
//        widomPump.setInterval(1L<<60);
        sim.integrator.getEventManager().addListener(widomPump);
        dataStreamPumps.add(widomPump);
        AccumulatorFactory muFactory = new AccumulatorFactory() {
            @Override
            public AccumulatorAverageFixed makeAccumulator() {
                return new AccumulatorMimicMu(sim.integrator);
            }
        };
        addAsTab(createStatPanel(widomFork, d, muFactory, false), "Chemical Potential", true);
        AccumulatorAverageCollapsing widomAvg = new AccumulatorAverageCollapsing();
        widomFork.addDataSink(widomAvg);
        AccumulatorMimicMu accMu = new AccumulatorMimicMu(sim.integrator);
        widomAvg.addDataSink(accMu);

        final DisplayTextBoxesCAE dDisplay = new DisplayTextBoxesCAE();
        dDisplay.setAccumulator(dAccumulator);
        final DisplayTextBoxesCAE pDisplay = new DisplayTextBoxesCAE();
        pDisplay.setAccumulator(pAccumulator);
        final DisplayTextBoxesCAE peDisplay = new DisplayTextBoxesCAE();
        peDisplay.setAccumulator(peAccumulator);
        final DisplayTextBoxesCAE muDisplay = new DisplayTextBoxesCAE();
        muDisplay.setAccumulator(accMu);
        muDisplay.setDoShowCurrent(false);

        //************* Lay out components ****************//

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        //temperature selector
        DeviceThermoSlider temperatureSelect = new DeviceThermoSlider(sim.getController(), sim.integrator);
        temperatureSelect.setPrecision(1);
        temperatureSelect.setMinimum(0.0);
        temperatureSelect.setMaximum(10.0);
        temperatureSelect.setSliderMajorValues(4);
        temperatureSelect.setIsothermalButtonsVisibility(false);

        IAction resetAction = new IAction() {
            public void actionPerformed() {
                meterWidom.setTemperature(sim.integrator.getTemperature());

                double dMin = ((int) (densityMeter.getDataAsScalar() * 100)) * 0.01;
                if (dMin > 0) dMin -= 0.005;

                // Reset density (Density is set and won't change, but
                // do this anyway)
                dDisplay.putData(dAccumulator.getData());
                dDisplay.repaint();

                // IS THIS WORKING?
                pPump.actionPerformed();
                pDisplay.putData(pAccumulator.getData());
                pDisplay.repaint();
                pePump.actionPerformed();
                peDisplay.putData(peAccumulator.getData());
                peDisplay.repaint();

                getDisplayBox(sim.box).graphic().repaint();
            }
        };

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

        add(dDisplay);
        add(pDisplay);
        add(peDisplay);
        add(muDisplay);

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

    protected DataFork makeHistoryPlot(ArrayList<DataPump> dataStreamPumps, DataSourceCountSteps timeCounter, JPanel historyPanel, AccumulatorAverageCollapsing pAccumulator, DataPumpListener pPump, String name) {
        AccumulatorHistory pHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        pHistory.setTimeDataSource(timeCounter);
        DataFork pFork = new DataFork(new IDataSink[]{pHistory, pAccumulator});
        pPump.setDataSink(pFork);
        sim.integrator.getEventManager().addListener(pPump);
        pAccumulator.setPushInterval(1);
        dataStreamPumps.add(pPump);

        DisplayPlot pPlot = new DisplayPlot();
        pHistory.setDataSink(pPlot.getDataSet().makeDataSink());
        pPlot.setDoLegend(false);
        pPlot.getPlot().setYLabel(name);
        pPlot.setLabel(name);
        historyPanel.add(pPlot.graphic());
        AccumulatorHistory[] pSinks = new AccumulatorHistory[3];
        for (int i = 0; i < pSinks.length; i++) {
            pSinks[i] = new AccumulatorHistory(new HistoryCollapsingDiscard());
            pSinks[i].setTimeDataSource(timeCounter);
        }
        DataProcessorBounds dpBounds = new DataProcessorBounds(pSinks);
        pAccumulator.addDataSink(dpBounds, new AccumulatorAverage.StatType[]{pAccumulator.AVERAGE, pAccumulator.ERROR});
        String[] labels = new String[]{"avg-", "avg", "avg+"};
        for (int i = 0; i < pSinks.length; i++) {
            pSinks[i].addDataSink(pPlot.getDataSet().makeDataSink());
            pPlot.setLegend(new DataTag[]{pSinks[i].getTag()}, labels[i]);
        }
        return pFork;
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
            if (doHistory && i % 5 == 0) {
                AccumulatorHistory accBlockHistory = new AccumulatorHistory(new HistoryStatistics(100));
                bit.setBlockDataSink(accBlockHistory);
                accBlockHistory.addDataSink(blockHistoryPlot.getDataSet().makeDataSink());
                accBlockHistory.setPushInterval(1 + 2000 / (1L << i));
                blockHistoryPlot.setLegend(new DataTag[]{accBlockHistory.getTag()}, "" + (1L << i));
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

    public static void main(String[] args) {
        Space sp = null;
        if (args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    sp = Space3D.getInstance();
                } else {
                    sp = Space2D.getInstance();
                }
            } catch (NumberFormatException e) {
            }
        } else {
            sp = Space3D.getInstance();
        }

        StatisticsMCGraphic ljmcGraphic = new StatisticsMCGraphic(new LJMC(sp), sp);
        SimulationGraphic.makeAndDisplayFrame
                (ljmcGraphic.getPanel(), APP_NAME);
    }

}


