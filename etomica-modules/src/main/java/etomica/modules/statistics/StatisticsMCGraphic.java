/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.statistics;

import etomica.action.IAction;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.graphics.*;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.modifier.ModifierBoolean;
import etomica.modules.ensembles.LJMC;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Quantity;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.List;

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
        dPlot.setLabel("Density");

        MeterPotentialEnergyFromIntegrator peMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorHistory peHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        peHistory.setTimeDataSource(timeCounter);
        final AccumulatorAverageCollapsing peAccumulator = new AccumulatorAverageCollapsing();
        peAccumulator.setPushInterval(10);
        DataFork peFork = new DataFork(new IDataSink[]{peHistory, peAccumulator});
        final DataPumpListener pePump = new DataPumpListener(peMeter, peFork, 1);
        sim.integrator.getEventManager().addListener(pePump);
        peHistory.setPushInterval(1);
        dataStreamPumps.add(pePump);

        DataCollector collectorErr = new DataCollector();
        DataCollector collectorCor = new DataCollector();
        for (int i = 0; i < 30; i++) {
            AccumulatorAverageFixed acc = new AccumulatorAverageFixed(1L << i);
            peFork.addDataSink(acc);
            acc.addDataSink(new DataAccPusher(i, collectorErr), new AccumulatorAverage.StatType[]{acc.ERROR});
            acc.addDataSink(new DataAccPusher(i, collectorCor), new AccumulatorAverage.StatType[]{acc.BLOCK_CORRELATION});
        }
        DisplayPlot peErrorPlot = new DisplayPlot();
        DataPumpListener peAccPumpErr = new DataPumpListener(collectorErr, peErrorPlot.getDataSet().makeDataSink(), 100);
        sim.integrator.getEventManager().addListener(peAccPumpErr);
        peErrorPlot.setLabel("PE Error");
        peErrorPlot.setDoLegend(false);
        peErrorPlot.getPlot().setXLog(true);
        add(peErrorPlot);
        DisplayPlot peCorPlot = new DisplayPlot();
        DataPumpListener peAccPumpCor = new DataPumpListener(collectorCor, peCorPlot.getDataSet().makeDataSink(), 100);
        sim.integrator.getEventManager().addListener(peAccPumpCor);
        peCorPlot.setLabel("PE Correlation");
        peCorPlot.setDoLegend(false);
        peCorPlot.getPlot().setXLog(true);
        peCorPlot.getPlot().setYRange(-1, 1);
        add(peCorPlot);

        DisplayPlot ePlot = new DisplayPlot();
        peHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setDoLegend(false);
        ePlot.setLabel("Energy");

        MeterPressure pMeter = new MeterPressure(space);
        pMeter.setIntegrator(sim.integrator);
        pMeter.setBox(sim.box);
        final AccumulatorAverageCollapsing pAccumulator = new AccumulatorAverageCollapsing();
        AccumulatorHistory pHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        pHistory.setTimeDataSource(timeCounter);
        DataFork pFork = new DataFork(new IDataSink[]{pHistory, pAccumulator});
        final DataPumpListener pPump = new DataPumpListener(pMeter, pFork, 1000);
        sim.integrator.getEventManager().addListener(pPump);
        pAccumulator.setPushInterval(1);
        dataStreamPumps.add(pPump);

        DisplayPlot pPlot = new DisplayPlot();
        pHistory.setDataSink(pPlot.getDataSet().makeDataSink());
        pPlot.setDoLegend(false);

        pPlot.setLabel("Pressure");


        final DisplayTextBoxesCAE dDisplay = new DisplayTextBoxesCAE();
        dDisplay.setAccumulator(dAccumulator);
        final DisplayTextBoxesCAE pDisplay = new DisplayTextBoxesCAE();
        pDisplay.setAccumulator(pAccumulator);
        final DisplayTextBoxesCAE peDisplay = new DisplayTextBoxesCAE();
        peDisplay.setAccumulator(peAccumulator);

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

        add(dPlot);
        add(ePlot);
        add(pPlot);
        add(dDisplay);
        add(pDisplay);
        add(peDisplay);

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

        StatisticsMCGraphic ljmdGraphic = new StatisticsMCGraphic(new LJMC(sp), sp);
        SimulationGraphic.makeAndDisplayFrame
                (ljmdGraphic.getPanel(), APP_NAME);
    }

    public static class Applet extends JApplet {

        public void init() {
            getRootPane().putClientProperty(
                    "defeatSystemEventQueueCheck", Boolean.TRUE);
            Space sp = Space3D.getInstance();
            StatisticsMCGraphic ljmdGraphic = new StatisticsMCGraphic(new LJMC(sp), sp);

            getContentPane().add(ljmdGraphic.getPanel());
        }

        private static final long serialVersionUID = 1L;
    }

    public static class AccumulatorAdder implements IntegratorListener {

        protected final DataFork fork;
        protected final List<AccumulatorAverageFixed> accumulators;
        protected final DataCollector collector;
        protected final AccumulatorAverageFixed.StatType stat;

        public AccumulatorAdder(DataFork fork, DataCollector collector) {
            this.fork = fork;
            accumulators = new ArrayList<>();
            this.collector = collector;
            this.stat = AccumulatorAverage.ERROR;
        }

        @Override
        public void integratorInitialized(IntegratorEvent e) {
            if (accumulators.size() == 0) {
                AccumulatorAverageFixed acc1 = new AccumulatorAverageFixed(1);
                fork.addDataSink(acc1);
                accumulators.add(acc1);
                acc1.addDataSink(new DataAccPusher(0, collector), new AccumulatorAverage.StatType[]{stat});
            }
        }

        @Override
        public void integratorStepStarted(IntegratorEvent e) {
        }

        @Override
        public void integratorStepFinished(IntegratorEvent e) {
            long n = accumulators.get(0).getSampleCount();
            int nAcc = 64 - Long.numberOfLeadingZeros(n);
            System.out.println(n + " " + nAcc);
            if (nAcc > accumulators.size()) {
                System.out.println("adding acc for bs=" + (1L << (nAcc - 1)));
                AccumulatorAverageFixed acc = new AccumulatorAverageFixed(1L << (nAcc - 1));
                fork.addDataSink(acc);
                accumulators.add(acc);
                acc.addDataSink(new DataAccPusher(accumulators.size() - 1, collector), new AccumulatorAverage.StatType[]{stat});
            }
        }
    }

    public static class DataAccPusher implements IDataSink {

        protected final int idx;
        protected final DataCollector collector;

        public DataAccPusher(int idx, DataCollector collector) {
            this.idx = idx;
            this.collector = collector;
        }

        @Override
        public void putData(IData data) {
            collector.setData(idx, data.getValue(0));
        }

        @Override
        public void putDataInfo(IEtomicaDataInfo dataInfo) {
            // don't care!
        }

        @Override
        public DataPipe getDataCaster(IEtomicaDataInfo inputDataInfo) {
            return null;
        }
    }

    public static class DataCollector implements IEtomicaDataSource {

        protected DataFunction data;
        protected DataFunction.DataInfoFunction dataInfo;
        protected final DataTag tag;

        public DataCollector() {
            tag = new DataTag();
            setLength(0);
        }

        public void setData(int i, double x) {
            if (data.getLength() <= i) {
                setLength(i + 1);
            }
            data.getData()[i] = x;
        }

        protected void setLength(int newLength) {
            double[] oldY = null;
            int oldSize = 0;
            if (dataInfo != null) {
                oldY = data.getData();
            }
            data = new DataFunction(new int[]{newLength});
            if (oldY != null) System.arraycopy(oldY, 0, data.getData(), 0, oldY.length);
            double[] xData = new double[newLength];
            for (int j = 0; j < newLength; j++) {
                xData[j] = 1L << j;
            }
            DataDoubleArray.DataInfoDoubleArray xDataInfo = new DataDoubleArray.DataInfoDoubleArray("block size", Quantity.DIMENSION, new int[]{newLength});
            dataInfo = new DataFunction.DataInfoFunction("stuff", Null.DIMENSION, new DataSourceIndependentSimple(xData, xDataInfo));
            dataInfo.addTag(tag);
        }

        @Override
        public DataTag getTag() {
            return tag;
        }

        @Override
        public IEtomicaDataInfo getDataInfo() {
            return dataInfo;
        }

        @Override
        public IData getData() {
            return data;
        }
    }
}


