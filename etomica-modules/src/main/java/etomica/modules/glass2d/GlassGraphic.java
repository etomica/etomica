/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.glass2d;

import etomica.action.IAction;
import etomica.atom.DiameterHashByType;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.*;
import etomica.data.types.DataDouble;
import etomica.data.types.DataTensor;
import etomica.graphics.*;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierBoolean;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Null;
import etomica.util.ParseArgs;

import java.awt.*;
import java.util.ArrayList;

public class GlassGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Lennard-Jones Molecular Dynamics";
    private final static int REPAINT_INTERVAL = 20;
    protected SimGlass sim;

    public GlassGraphic(final SimGlass simulation) {

        super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

        this.sim = simulation;

        sim.activityIntegrate.setSleepPeriod(1);

        //display of box, timer
        ColorSchemeByType colorScheme = new ColorSchemeByType();
        colorScheme.setColor(sim.speciesA.getLeafType(), Color.red);
        colorScheme.setColor(sim.speciesB.getLeafType(), Color.blue);
//        DisplayBox dbox0 = getDisplayBox(sim.box);
//        dbox0.setColorScheme(colorScheme);
        DiameterHashByType diameterHash = (DiameterHashByType) getDisplayBox(sim.box).getDiameterHash();
        diameterHash.setDiameter(sim.speciesB.getLeafType(), sim.potentialChoice == SimGlass.PotentialChoice.LJ ? 0.88 : 1 / 1.4);
        diameterHash.setDiameter(sim.speciesA.getLeafType(), 1);
//        sim.integrator.addListener(new IntervalActionAdapter(this.getDisplayBoxPaintAction(sim.box)));

        //meters and displays
        final MeterRDF rdfMeter = new MeterRDF(sim.getSpace());
        IntegratorListenerAction rdfMeterListener = new IntegratorListenerAction(rdfMeter);
        sim.integrator.getEventManager().addListener(rdfMeterListener);
        rdfMeterListener.setInterval(10);
        rdfMeter.getXDataSource().setXMax(4.0);
        rdfMeter.setBox(sim.box);
        DisplayPlot rdfPlot = new DisplayPlot();
        DataPump rdfPump = new DataPump(rdfMeter, rdfPlot.getDataSet().makeDataSink());
        IntegratorListenerAction rdfPumpListener = new IntegratorListenerAction(rdfPump);
        sim.integrator.getEventManager().addListener(rdfPumpListener);
        rdfPumpListener.setInterval(10);
        dataStreamPumps.add(rdfPump);

        rdfPlot.setDoLegend(false);
        rdfPlot.getPlot().setTitle("Radial Distribution Function");
        rdfPlot.setLabel("RDF");

        DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);

        ConfigurationStorage configStorageMSD = new ConfigurationStorage(sim.box, true);
        configStorageMSD.setEnabled(false); // start isothermal
        sim.integrator.getEventManager().addListener(configStorageMSD);

        ConfigurationStorage configStorage = new ConfigurationStorage(sim.box, false);
        configStorage.setEnabled(false); // start isothermal
        sim.integrator.getEventManager().addListener(configStorage);
        DisplayBox dbox = new DisplayBox(sim, sim.box);
        DisplayBoxCanvas2DGlass canvas = new DisplayBoxCanvas2DGlass(dbox, sim.getSpace(), sim.getController(), configStorage);
        remove(getDisplayBox(sim.box));
        dbox.setBoxCanvas(canvas);
        add(dbox);
        dbox.setColorScheme(colorScheme);
        dbox.setDiameterHash(diameterHash);
        canvas.setVisible(false);
        canvas.setVisible(true);

        ColorSchemeDeviation colorSchemeDeviation = new ColorSchemeDeviation(sim.box, configStorage);

        DataSourcePrevTime dsPrevTime = new DataSourcePrevTime(configStorage);
        DeviceSlider prevConfigSlider = new DeviceSlider(sim.getController(), new Modifier() {
            @Override
            public void setValue(double newValue) {
                int idx = (int) Math.round(newValue);
                canvas.setConfigIndex(idx);
                colorSchemeDeviation.setConfigIndex(idx);
                dsPrevTime.setPrevConfigIndex(idx);
            }

            @Override
            public double getValue() {
                return canvas.getConfigIndex();
            }

            @Override
            public Dimension getDimension() {
                return Null.DIMENSION;
            }

            @Override
            public String getLabel() {
                return "previous config (log2)";
            }
        });
        prevConfigSlider.setMaximum(30);
        prevConfigSlider.setNMajor(5);
        prevConfigSlider.setShowBorder(true);
        add(prevConfigSlider);
        DisplayTextBox displayPrevTime = new DisplayTextBox();
        DataPumpListener pumpPrevTime = new DataPumpListener(dsPrevTime, displayPrevTime, 1);
        sim.integrator.getEventManager().addListener(pumpPrevTime);
        add(displayPrevTime);

        DeviceCheckBox colorCheckbox = new DeviceCheckBox("color by displacement", new ModifierBoolean() {
            @Override
            public void setBoolean(boolean b) {
                if ((dbox.getColorScheme() == colorSchemeDeviation) == b) return;
                if (b) {
                    if (sim.integrator.isIsothermal()) return;
                    dbox.setColorScheme(colorSchemeDeviation);
                } else {
                    dbox.setColorScheme(colorScheme);
                }
            }

            @Override
            public boolean getBoolean() {
                return dbox.getColorScheme() == colorSchemeDeviation;
            }
        });
        colorCheckbox.setController(sim.getController());
        add(colorCheckbox);

        DeviceCheckBox showDispCheckbox = new DeviceCheckBox("show displacement", new ModifierBoolean() {
            @Override
            public void setBoolean(boolean b) {
                if (canvas.getDrawDisplacement() == b || sim.integrator.isIsothermal()) return;
                canvas.setDrawDisplacement(b);
            }

            @Override
            public boolean getBoolean() {
                return canvas.getDrawDisplacement();
            }
        });

        showDispCheckbox.setController(sim.getController());
        add(showDispCheckbox);




        IAction repaintAction = new IAction() {
            public void actionPerformed() {
                dbox.repaint();
            }
        };
        IntegratorListenerAction repaintAction2 = new IntegratorListenerAction(repaintAction);
        repaintAction2.setInterval(100);
        sim.integrator.getEventManager().addListener(repaintAction2);

        final IAction resetDataAction = new IAction() {
            public void actionPerformed() {
                getController().getSimRestart().getDataResetAction().actionPerformed();
                rdfMeter.reset();
            }
        };


        //add meter and display for current kinetic temperature

        MeterTemperature thermometer = new MeterTemperature(sim.box, space.D());
        DataFork temperatureFork = new DataFork();
        final DataPump temperaturePump = new DataPump(thermometer, temperatureFork);
        IntegratorListenerAction temperaturePumpListener = new IntegratorListenerAction(temperaturePump);
        sim.integrator.getEventManager().addListener(temperaturePumpListener);
        temperaturePumpListener.setInterval(10);
        final AccumulatorHistory temperatureHistory = new AccumulatorHistory();
        temperatureHistory.setTimeDataSource(timeCounter);
        temperatureFork.setDataSinks(new IDataSink[]{temperatureHistory});

        dataStreamPumps.add(temperaturePump);

        // Number density box
        MeterDensity densityMeter = new MeterDensity(sim.getSpace());
        densityMeter.setBox(sim.box);
        final DisplayTextBox densityBox = new DisplayTextBox();
        final DataPump densityPump = new DataPump(densityMeter, densityBox);
        IntegratorListenerAction densityPumpListener = new IntegratorListenerAction(densityPump);
        sim.integrator.getEventManager().addListener(densityPumpListener);
        densityPumpListener.setInterval(10);
        dataStreamPumps.add(densityPump);
        densityBox.setLabel("Number Density");

        AccumulatorHistory energyHistory = null, peHistory = null;
        final AccumulatorAverageCollapsing peAccumulator;
        DataFork peFork = null;
        if (sim.potentialChoice != SimGlass.PotentialChoice.HS) {
            MeterEnergy eMeter = new MeterEnergy(sim.integrator.getPotentialMaster(), sim.box);
            energyHistory = new AccumulatorHistory();
            energyHistory.setTimeDataSource(timeCounter);
            DataPump energyPump = new DataPump(eMeter, energyHistory);
            IntegratorListenerAction energyPumpListener = new IntegratorListenerAction(energyPump);
            sim.integrator.getEventManager().addListener(energyPumpListener);
            energyPumpListener.setInterval(60);
            energyHistory.setPushInterval(5);
            dataStreamPumps.add(energyPump);

            MeterPotentialEnergy peMeter = new MeterPotentialEnergy(sim.integrator.getPotentialMaster(), sim.box);
            peHistory = new AccumulatorHistory();
            peHistory.setTimeDataSource(timeCounter);
            peAccumulator = new AccumulatorAverageCollapsing();
            peAccumulator.setPushInterval(10);
            peFork = new DataFork(new IDataSink[]{peHistory, peAccumulator});
            DataPump pePump = new DataPump(peMeter, peFork);
            IntegratorListenerAction pePumpListener = new IntegratorListenerAction(pePump);
            sim.integrator.getEventManager().addListener(pePumpListener);
            pePumpListener.setInterval(60);
            peHistory.setPushInterval(5);
            dataStreamPumps.add(pePump);
        } else {
            peAccumulator = null;
        }

        MeterKineticEnergy keMeter = new MeterKineticEnergy(sim.box);
        AccumulatorHistory keHistory = new AccumulatorHistory();
        keHistory.setTimeDataSource(timeCounter);
        DataFork keFork = new DataFork();
        DataPump kePump = new DataPump(keMeter, keFork);
        keFork.addDataSink(keHistory);
        final AccumulatorAverage keAvg = new AccumulatorAverageCollapsing();
        keFork.addDataSink(keAvg);
        IntegratorListenerAction kePumpListener = new IntegratorListenerAction(kePump);
        sim.integrator.getEventManager().addListener(kePumpListener);
        kePumpListener.setInterval(60);
        keHistory.setPushInterval(5);
        dataStreamPumps.add(kePump);

        DisplayPlot ePlot = new DisplayPlot();
        if (sim.potentialChoice != SimGlass.PotentialChoice.HS) {
            energyHistory.setDataSink(ePlot.getDataSet().makeDataSink());
            ePlot.setLegend(new DataTag[]{energyHistory.getTag()}, "Total");
            peHistory.setDataSink(ePlot.getDataSet().makeDataSink());
            ePlot.setLegend(new DataTag[]{peHistory.getTag()}, "Potential");
        }
        keHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{keHistory.getTag()}, "Kinetic");

        ePlot.getPlot().setTitle("Energy History");
        ePlot.setDoLegend(true);
        ePlot.setLabel("Energy");

        IDataSource pMeter;
        if (sim.integrator instanceof IntegratorVelocityVerlet) {
            pMeter = new MeterPressureTensorFromIntegrator(space);
            ((MeterPressureTensorFromIntegrator) pMeter).setIntegrator((IntegratorVelocityVerlet) sim.integrator);
        } else {
            pMeter = new MeterPressureHardTensor(space);
            ((MeterPressureHardTensor) pMeter).setIntegrator((IntegratorHard) sim.integrator);
            new MeterPressureHard((IntegratorHard) sim.integrator);
        }
        DataFork pTensorFork = new DataFork();
        DataPumpListener pPump = new DataPumpListener(pMeter, pTensorFork);
        AccumulatorAutocorrelationPTensor dpAutocor = new AccumulatorAutocorrelationPTensor(100000, sim.integrator.getTimeStep());
        pTensorFork.addDataSink(dpAutocor);
        dpAutocor.setPushInterval(1000);
        DataProcessorTensorTrace tracer = new DataProcessorTensorTrace();
        pTensorFork.addDataSink(tracer);
        DataFork pFork = new DataFork();
        tracer.setDataSink(pFork);
        sim.integrator.getEventManager().addListener(pPump);
        final AccumulatorAverageCollapsing pAccumulator = new AccumulatorAverageCollapsing();
        pFork.addDataSink(pAccumulator);
        pAccumulator.setPushInterval(100);
        dataStreamPumps.add(pPump);
        AccumulatorHistory historyP = new AccumulatorHistory(new HistoryCollapsingAverage());
        historyP.setTimeDataSource(timeCounter);
        pFork.addDataSink(historyP);
        DisplayPlot plotP = new DisplayPlot();
        historyP.addDataSink(plotP.getDataSet().makeDataSink());
        plotP.setLabel("P");
        add(plotP);

        DisplayPlot plotPTensorAutocor = new DisplayPlot();
        plotPTensorAutocor.setLabel("P Tensor autocor");
        plotPTensorAutocor.setLegend(new DataTag[]{dpAutocor.getTag()}, "avg");
        dpAutocor.addDataSink(plotPTensorAutocor.getDataSet().makeDataSink());
        add(plotPTensorAutocor);

        DisplayTextBoxesCAE peDisplay = null;
        if (sim.potentialChoice != SimGlass.PotentialChoice.HS) {
            AccumulatorHistory historyPE = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyPE.setTimeDataSource(timeCounter);
            peFork.addDataSink(historyPE);
            DisplayPlot plotPE = new DisplayPlot();
            historyPE.addDataSink(plotPE.getDataSet().makeDataSink());
            plotPE.setLabel("PE");
            add(plotPE);
            peDisplay = new DisplayTextBoxesCAE();
            peDisplay.setAccumulator(peAccumulator);
        }

        final DisplayTextBoxesCAE pDisplay = new DisplayTextBoxesCAE();
        pDisplay.setAccumulator(pAccumulator);

        DataSourceMSD meterMSD = new DataSourceMSD(configStorageMSD);
        configStorageMSD.addListener(meterMSD);
        DisplayPlot plotMSD = new DisplayPlot();
        DataPumpListener pumpMSD = new DataPumpListener(meterMSD, plotMSD.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpMSD);
        plotMSD.setLabel("MSD");
        plotMSD.getPlot().setYLog(true);
        plotMSD.getPlot().setXLog(true);
        add(plotMSD);

        //************* Lay out components ****************//

        DeviceCheckBox swapCheckbox = new DeviceCheckBox("isothermal", new ModifierBoolean() {
            @Override
            public void setBoolean(boolean b) {

                if (sim.integrator.isIsothermal() == b) return;
                if (b) {
                    sim.integrator.setIsothermal(true);
                    sim.integrator.setIntegratorMC(sim.integratorMC, 10000);
                    dbox.setColorScheme(colorScheme);
                    canvas.setDrawDisplacement(false);
                    configStorage.reset();
                    configStorage.setEnabled(false);
                    configStorageMSD.reset();
                    configStorageMSD.setEnabled(false);
                    dpAutocor.reset();
                } else {
                    sim.integrator.setIntegratorMC(null, 0);
                    sim.integrator.setIsothermal(false);
                    configStorage.reset();
                    configStorage.setEnabled(true);
                    configStorageMSD.reset();
                    configStorageMSD.setEnabled(true);
                    if (colorCheckbox.getState()) dbox.setColorScheme(colorSchemeDeviation);
                    if (showDispCheckbox.getState()) canvas.setDrawDisplacement(true);
                    dpAutocor.reset();
                }
            }

            @Override
            public boolean getBoolean() {
                return sim.integrator.isIsothermal();
            }
        });
        swapCheckbox.setController(sim.getController());
        add(swapCheckbox);

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        getDisplayBox(sim.box).setScale(0.7);

        //temperature selector
        DeviceSlider temperatureSelect = new DeviceSlider(sim.getController(), sim.integrator, "temperature");
        temperatureSelect.setPrecision(2);
        temperatureSelect.setNMajor(4);
        temperatureSelect.setMinimum(0.0);
        temperatureSelect.setMaximum(2.0);
        temperatureSelect.setLabel("Temperature");
        temperatureSelect.setShowBorder(true);
        temperatureSelect.setShowValues(true);
        temperatureSelect.setEditValues(true);

        temperatureSelect.setPostAction(new IAction() {
            public void actionPerformed() {
                resetDataAction.actionPerformed();
                getDisplayBox(sim.box).repaint();
            }
        });
        configStorage.reset();
        dbox.repaint();

        IAction resetAction = new IAction() {
            public void actionPerformed() {
                rdfMeter.reset();

                getDisplayBox(sim.box).graphic().repaint();
            }
        };

        this.getController().getReinitButton().setPostAction(resetAction);
        this.getController().getResetAveragesButton().setPostAction(resetAction);

        getPanel().controlPanel.add(temperatureSelect.graphic(), vertGBC);

        add(rdfPlot);
        add(ePlot);
        add(densityBox);
        add(pDisplay);
        if (sim.potentialChoice != SimGlass.PotentialChoice.HS) {
            add(peDisplay);
        }

    }

    public static void main(String[] args) {
        SimGlass.GlassParams params = new SimGlass.GlassParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.doSwap = true;
            params.potential = SimGlass.PotentialChoice.SS;
            params.nA = params.nB = 400;
            params.density = 1.35;
        }
        SimGlass sim = new SimGlass(params.D, params.nA, params.nB, params.density, params.doSwap, params.potential);

        GlassGraphic ljmdGraphic = new GlassGraphic(sim);
        SimulationGraphic.makeAndDisplayFrame
                (ljmdGraphic.getPanel(), APP_NAME);
    }

    /**
     * Inner class to find the total pressure of the system from the pressure
     * tensor.
     */
    public static class DataProcessorTensorTrace extends DataProcessor {

        public DataProcessorTensorTrace() {
            data = new DataDouble();
        }

        protected IData processData(IData inputData) {
            // take the trace and divide by the dimensionality
            data.x = ((DataTensor) inputData).x.trace() / ((DataTensor) inputData).x.D();
            return data;
        }

        protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
            if (!(inputDataInfo instanceof DataTensor.DataInfoTensor)) {
                throw new IllegalArgumentException("Gotta be a DataInfoTensor");
            }
            dataInfo = new DataDouble.DataInfoDouble(inputDataInfo.getLabel(), inputDataInfo.getDimension());
            return dataInfo;
        }

        protected final DataDouble data;
    }

}


