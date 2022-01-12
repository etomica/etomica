/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.ensembles;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.DiameterHash;
import etomica.atom.DiameterHashByType;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.graphics.*;
import etomica.normalmode.MeterPressureHMA;
import etomica.space.Space;
import etomica.space1d.Space1D;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;

import java.awt.*;
import java.io.IOException;
import java.util.ArrayList;

public class LJMCGraphic1D extends SimulationGraphic {

    private final static String APP_NAME = "Ensembles";
    private final static int REPAINT_INTERVAL = 100;
    protected final LJMC1D sim;

    public LJMCGraphic1D(final LJMC1D simulation) {
        // simulation - this contains a LJMC1D object. This object is created once.
        // _space - this variable is a space object.

        // Constructed from the base-class SimulationGraphic
        super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);

        // An arraylist containing DataPump objects.
        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

        this.sim = simulation;      // creating an instance of the simulation.

        // Display of box, timer
        ColorSchemeByType colorScheme = new ColorSchemeByType();
        colorScheme.setColor(sim.species.getLeafType(), Color.red);
        getDisplayBox(sim.box).setColorScheme(new ColorSchemeByType());
        ((DiameterHashByType)getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species.getLeafType(), 0.5);

        DataSourceCountSteps timeCounter = new DataSourceCountSteps(sim.integrator);

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
        final MeterPressureHMA pMeterHMA;
        try {
            pMeterHMA = new MeterPressureHMA(space, sim.potentialMaster, sim.coordinates, true, 0);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
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

        add(ePlot);
        add(pPlot);
        add(pDisplay);
        add(pDisplayHMA);
        add(pDisplayConventional);
        add(peDisplay);

        final DeviceButton slowButton = new DeviceButton(sim.getController(), null);
        slowButton.setAction(new IAction() {
            public void actionPerformed() {
                int sleep = (int) sim.getController().getSleepPeriod();
                sleep = 1-sleep;
                sim.getController().setSleepPeriod(sleep);
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
                } else {
                    sp = Space2D.getInstance();
                }
            } catch (NumberFormatException e) {
            }
        } else {
            sp = Space1D.getInstance();
        }

        LJMC1D sim = new LJMC1D(sp, 10.0, 20, 3.0, 1.5, new int[]{8});
        sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
        LJMCGraphic1D ljmcGraphic = new LJMCGraphic1D(sim);
        SimulationGraphic.makeAndDisplayFrame
                (ljmcGraphic.getPanel(), APP_NAME);
    }
}