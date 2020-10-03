/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.ensembles;

import etomica.action.IAction;
import etomica.data.*;
import etomica.data.histogram.HistogramCollapsing;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.*;
import etomica.graphics.*;
import etomica.math.DoubleRange;
import etomica.modifier.ModifierBoolean;
import etomica.space.Space;
import etomica.space1d.Space1D;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.normalmode.MeterPressureHMA;
import etomica.data.DataSplitter;

import javax.swing.*;
import java.awt.*;
import java.io.IOException;
import java.util.ArrayList;

public class LJMCGraphic1D extends SimulationGraphic {

    private final static String APP_NAME = "Ensembles";
    private final static int REPAINT_INTERVAL = 100;
    protected final LJMC1D sim;

    protected boolean volumeChanges = false;
    protected boolean constMu = false;

    public LJMCGraphic1D(final LJMC1D simulation, Space _space) throws IOException {
        // simulation - this contains a LJMC1D object. This object is created once.
        // _space - this variable is a space object.

        // Constructed from the base-class SimulationGraphic
        super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);

        // An arraylist containing DataPump objects.
        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

        this.sim = simulation;      // creating an instance of the simulation.

        // Display of box, timer
        ColorSchemeByType colorScheme = new ColorSchemeByType();
        colorScheme.setColor(sim.species.getLeafType(),java.awt.Color.red);
        getDisplayBox(sim.box).setColorScheme(new ColorSchemeByType());

        DataSourceCountSteps timeCounter = new DataSourceCountSteps(sim.integrator);

        // Number density box
        final MeterDensity densityMeter = new MeterDensity(sim.getSpace());
        densityMeter.setBox(sim.box);
        AccumulatorHistory dHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        dHistory.setTimeDataSource(timeCounter);
        final AccumulatorAverageCollapsing dAccumulator = new AccumulatorAverageCollapsing();
        dAccumulator.setPushInterval(10);
        final HistogramCollapsing dh = new HistogramCollapsing();
        double dmin = ((int)(densityMeter.getDataAsScalar() * 100))*0.01;
        if (dmin > 0) dmin -= 0.005;
        dh.setXRange(new DoubleRange(dmin, dmin+0.1));
        final AccumulatorHistogram dHistogram = new AccumulatorHistogram(dh);
        dHistogram.setPushInterval(10);
        DataFork dFork = new DataFork(new IDataSink[]{dHistory, dAccumulator, dHistogram});
        final DataPumpListener dPump = new DataPumpListener(densityMeter, dFork, 100);
        sim.integrator.getEventManager().addListener(dPump);
        dHistory.setPushInterval(1);
        dataStreamPumps.add(dPump);

        DisplayPlotXChart dPlot = new DisplayPlotXChart();
        dHistory.setDataSink(dPlot.getDataSet().makeDataSink());
        dPlot.setDoLegend(false);
        dPlot.setLabel("Density");

        DisplayPlotXChart dHistogramPlot = new DisplayPlotXChart();
        dHistogram.setDataSink(dHistogramPlot.getDataSet().makeDataSink());
        dHistogramPlot.setDoLegend(false);
        dHistogramPlot.setLabel("Density histogram");

        // Below is a block of code, corresponding to the Potential Energy.
        MeterPotentialEnergyFromIntegrator peMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorHistory peHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        peHistory.setTimeDataSource(timeCounter);
        final AccumulatorAverageCollapsing peAccumulator = new AccumulatorAverageCollapsing();
        peAccumulator.setPushInterval(10);
        DataFork peFork = new DataFork(new IDataSink[]{peHistory, peAccumulator});  // This is the creation of data sink array; two sinks - peHistory and peAccumulator.
        final DataPumpListener pePump = new DataPumpListener(peMeter, peFork, 100); // Takes data from source (peMeter) to sink (peFork).
        sim.integrator.getEventManager().addListener(pePump);
        peHistory.setPushInterval(1);
        dataStreamPumps.add(pePump);

        // Below is a block of code that is responsible for displaying the Energy.
        DisplayPlotXChart ePlot = new DisplayPlotXChart();
        peHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setDoLegend(false);
        ePlot.setLabel("Energy");

        // Below is a block of code, corresponding to the  Conventional Pressure.
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

        // Below is a block of code corresponding to the Harmonically-Mapped Average Pressure.
        // TODO: Clean up this section. Does every single method need to be called right here?
        MeterPressureHMA pMeterHMA = new MeterPressureHMA(space, sim.potentialMaster, sim.coordinates, true);
        pMeterHMA.setTemperature(sim.integrator.getTemperature());
        pMeterHMA.setTruncationRadius(4.1); // Avoid using multiples of the lattice constant.
        pMeterHMA.setPRes();
        pMeterHMA.calculateGij();

        AccumulatorHistory pHistoryHMA = new AccumulatorHistory(new HistoryCollapsingAverage());
        pHistoryHMA.setTimeDataSource(timeCounter);
        AccumulatorHistory pHistoryHMANormalMode = new AccumulatorHistory(new HistoryCollapsingAverage());
        pHistoryHMANormalMode.setTimeDataSource(timeCounter);
        AccumulatorHistory pHistoryQuasi = new AccumulatorHistory(new HistoryCollapsingAverage());
        pHistoryQuasi.setTimeDataSource(timeCounter);
        AccumulatorHistory pHistoryConventional = new AccumulatorHistory(new HistoryCollapsingAverage());
        pHistoryConventional.setTimeDataSource(timeCounter);


        final AccumulatorAverageCollapsing pAccumulatorHMARealSpace = new AccumulatorAverageCollapsing();
        final AccumulatorAverageCollapsing pAccumulatorQuasi = new AccumulatorAverageCollapsing();
        final AccumulatorAverageCollapsing pAccumulatorHMANormalMode = new AccumulatorAverageCollapsing();
        final AccumulatorAverageCollapsing pAccumulatorConventional = new AccumulatorAverageCollapsing();

        DataSplitter splitter = new DataSplitter();
        final DataPumpListener pPumpHMAQuasi = new DataPumpListener(pMeterHMA, splitter, 1000);
        DataFork pForkQuasi = new DataFork(new IDataSink[]{pHistoryQuasi, pAccumulatorQuasi});
        DataFork pForkHMA = new DataFork(new IDataSink[]{pHistoryHMA, pAccumulatorHMARealSpace});
        DataFork pForkHMANormalMode = new DataFork(new IDataSink[]{pHistoryHMANormalMode, pAccumulatorHMANormalMode});
        DataFork pForkConventional = new DataFork(new IDataSink[]{pHistoryConventional, pAccumulatorConventional});
        splitter.setDataSink(1, pForkConventional);
        splitter.setDataSink(3, pForkQuasi);
        splitter.setDataSink(4, pForkHMANormalMode);
        splitter.setDataSink(5, pForkHMA);

        sim.integrator.getEventManager().addListener(pPumpHMAQuasi);
        pAccumulatorHMARealSpace.setPushInterval(1);
        pAccumulatorQuasi.setPushInterval(1);
        pAccumulatorHMANormalMode.setPushInterval(1);
        pAccumulatorConventional.setPushInterval(1);
        dataStreamPumps.add(pPumpHMAQuasi);

        // Below is a block of code that is responsible for displaying the Pressure.
        DisplayPlotXChart pPlot = new DisplayPlotXChart();      // Creating a Pressure Plot Object.
        pHistory.setDataSink(pPlot.getDataSet().makeDataSink());        // This is for Conventional Pressure.
        pHistoryHMA.setDataSink(pPlot.getDataSet().makeDataSink());     // This is for HMA Pressure from real-space coordinates.
        pHistoryQuasi.setDataSink(pPlot.getDataSet().makeDataSink());       // This is for the QuasiHarmonic Pressure.
        pHistoryHMANormalMode.setDataSink(pPlot.getDataSet().makeDataSink());   // This is for the HMA Pressure from normal-mode coordinates.
        pHistoryConventional.setDataSink(pPlot.getDataSet().makeDataSink());        // This is for the Conventional Pressure from MeterHMA object.
        pPlot.setDoLegend(true);
        pPlot.setLegend(new DataTag[]{pHistoryHMA.getTag()}, "HMA - Real Space Coordinates");
        pPlot.setLegend(new DataTag[]{pHistoryQuasi.getTag()}, "Quasi-harmonic Pressure");
        pPlot.setLegend(new DataTag[]{pHistory.getTag()}, "Conventional from MeterPressure");
        pPlot.setLegend(new DataTag[]{pHistoryHMANormalMode.getTag()}, "HMA - Normal Mode Coordinates");
        pPlot.setLegend(new DataTag[]{pHistoryConventional.getTag()}, "Conventional from MeterPressureHMA");

        pPlot.setLabel("Pressure");

        final DisplayTextBoxesCAE dDisplay = new DisplayTextBoxesCAE();
        dDisplay.setAccumulator(dAccumulator);
        final DisplayTextBoxesCAE pDisplay = new DisplayTextBoxesCAE();
        pDisplay.setAccumulator(pAccumulator);
        final DisplayTextBoxesCAE peDisplay = new DisplayTextBoxesCAE();
        peDisplay.setAccumulator(peAccumulator);
        final DisplayTextBoxesCAE pDisplayHMA= new DisplayTextBoxesCAE();
        pDisplayHMA.setAccumulator(pAccumulatorHMARealSpace);
        final DisplayTextBoxesCAE pDisplayQuasi = new DisplayTextBoxesCAE();
        pDisplayQuasi.setAccumulator(pAccumulatorQuasi);
        final DisplayTextBoxesCAE pDisplayNormalModeHMA = new DisplayTextBoxesCAE();
        pDisplayNormalModeHMA.setAccumulator(pAccumulatorHMANormalMode);
        final DisplayTextBoxesCAE pDisplayConventional = new DisplayTextBoxesCAE();
        pDisplayConventional.setAccumulator(pAccumulatorConventional);

        //************* Lay out components ****************//

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        //temperature selector
        double T = sim.integrator.getTemperature();
        DeviceThermoSlider temperatureSelect = new DeviceThermoSlider(sim.getController(), sim.integrator);
        temperatureSelect.setPostAction(new IAction() {
            @Override
            public void actionPerformed() {
                pMeterHMA.setTemperature(temperatureSelect.getTemperature());
            }
        });
        temperatureSelect.setPrecision(1);
        temperatureSelect.setMinimum(0.0);
        temperatureSelect.setMaximum(10.0);
        temperatureSelect.setSliderMajorValues(4);
        temperatureSelect.setIsothermalButtonsVisibility(false);
        temperatureSelect.setTemperature(T);

        IAction resetAction = new IAction() {
            public void actionPerformed() {

                double dMin = ((int)(densityMeter.getDataAsScalar() * 100))*0.01;
                if (dMin > 0) dMin -= 0.005;
                dh.setXRange(new DoubleRange(dMin, dMin+0.1));

                // Reset density (Density is set and won't change, but
                // do this anyway)
                dPump.actionPerformed();
                dDisplay.putData(dAccumulator.getData());

                // IS THIS WORKING?
                pPump.actionPerformed();
                pPumpHMAQuasi.actionPerformed();
                pDisplay.putData(pAccumulator.getData());
                pDisplayHMA.putData(pAccumulatorHMARealSpace.getData());
                pDisplayQuasi.putData(pAccumulatorQuasi.getData());
                pDisplayNormalModeHMA.putData(pAccumulatorHMANormalMode.getData());
                pDisplayConventional.putData(pAccumulatorConventional.getData());
                pePump.actionPerformed();
                peDisplay.putData(peAccumulator.getData());

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
                }
                else {
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
                }
                else {
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
        add(dHistogramPlot);
        add(ePlot);
        add(pPlot);
        add(dDisplay);
        add(pDisplay);
        add(pDisplayHMA);
        add(pDisplayConventional);
        add(peDisplay);

        final DeviceButton slowButton = new DeviceButton(sim.getController(), null);
        slowButton.setAction(new IAction() {
            public void actionPerformed() {
                int sleep = (int) sim.activityIntegrate.getSleepPeriod();
                sleep = 1-sleep;
                sim.activityIntegrate.setSleepPeriod(sleep);
                slowButton.setLabel(sleep == 0 ? "Slow" : "Fast");
            }
        });
        slowButton.setLabel("Slow");
        add(slowButton);
    }

    public static void main(String[] args) {
        // By default, the created space will be 1D.
        Space sp = null;
        if(args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    sp = Space3D.getInstance();
                }
                else {
                    sp = Space2D.getInstance();
                }
            } catch(NumberFormatException e) {}
        }
        else {
            sp = Space1D.getInstance();
        }

        try {
            LJMCGraphic1D ljmcGraphic = new LJMCGraphic1D(new LJMC1D(sp), sp);
            SimulationGraphic.makeAndDisplayFrame
                    (ljmcGraphic.getPanel(), APP_NAME);
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static class Applet extends javax.swing.JApplet {

        public void init() {
            getRootPane().putClientProperty(
                    "defeatSystemEventQueueCheck", Boolean.TRUE);
            Space sp = Space1D.getInstance();
            LJMCGraphic1D ljmcGraphic = null;
            try {
                ljmcGraphic = new LJMCGraphic1D(new LJMC1D(sp), sp);
            } catch (IOException e) {
                e.printStackTrace();
            }

            getContentPane().add(ljmcGraphic.getPanel());
        }

        private static final long serialVersionUID = 1L;
    }
}