/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.glass;

import etomica.action.IAction;
import etomica.cavity.DataProcessorErrorBar;
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
import etomica.units.SimpleUnit;
import etomica.units.Unit;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;
import etomica.util.ParseArgs;

import javax.swing.*;
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

        Unit peUnit = new SimpleUnit(Energy.DIMENSION, sim.box.getLeafList().size(), "u", "u", false);

        //display of box, timer
        ColorSchemeByType colorScheme = new ColorSchemeByType();
        colorScheme.setColor(sim.speciesA.getLeafType(), Color.red);
        colorScheme.setColor(sim.speciesB.getLeafType(), Color.blue);

        //meters and displays
        final MeterRDF rdfMeter = new MeterRDF(sim.getSpace());
        IntegratorListenerAction rdfMeterListener = new IntegratorListenerAction(rdfMeter);
        sim.integrator.getEventManager().addListener(rdfMeterListener);
        rdfMeterListener.setInterval(100);
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

        ConfigurationStorage configStorageLinear = new ConfigurationStorage(sim.box, ConfigurationStorage.StorageType.LINEAR, 1024, 32);
        configStorageLinear.setEnabled(false); // start isothermal
        sim.integrator.getEventManager().addListener(configStorageLinear);

        ConfigurationStorage configStorageMSD = new ConfigurationStorage(sim.box, ConfigurationStorage.StorageType.MSD);
        configStorageMSD.setEnabled(false); // start isothermal
        sim.integrator.getEventManager().addListener(configStorageMSD);

        ConfigurationStorage configStorage = new ConfigurationStorage(sim.box, ConfigurationStorage.StorageType.LOG2);
        configStorage.setEnabled(false); // start isothermal
        sim.integrator.getEventManager().addListener(configStorage);
        DisplayBox dbox;
        DisplayBoxCanvasGlass canvas;
        DisplayCanvas c;
        dbox = new DisplayBox(sim, sim.box);
        if (sim.getSpace().D() == 2) {
            c = new DisplayBoxCanvas2DGlass(dbox, sim.getSpace(), sim.getController(), configStorage);
        } else {
            c = new DisplayBoxCanvas3DGlass(dbox, sim.getSpace(), sim.getController(), configStorage);
        }
        canvas = (DisplayBoxCanvasGlass) c;
        remove(getDisplayBox(sim.box));
        dbox.setBoxCanvas(c);
        add(dbox);
        c.setVisible(false);
        c.setVisible(true);
        dbox.setColorScheme(colorScheme);

        DiameterHashGlass diameterHash = new DiameterHashGlass();
        diameterHash.setDiameter(sim.speciesB.getLeafType(), sim.potentialChoice == SimGlass.PotentialChoice.LJ ? 0.88 : 1 / 1.4);
        diameterHash.setDiameter(sim.speciesA.getLeafType(), 1);
        dbox.setDiameterHash(diameterHash);

        AtomTestDeviation atomFilterDeviation = new AtomTestDeviation(sim.box, configStorage);

        ColorSchemeDeviation colorSchemeDeviation = new ColorSchemeDeviation(sim.box, configStorage);
        ColorSchemeDirection colorSchemeDirection = new ColorSchemeDirection(sim.box, configStorage);

        DataSourcePrevTime dsPrevTime = new DataSourcePrevTime(configStorage);


        DeviceSlider prevConfigSlider = new DeviceSlider(sim.getController(), new Modifier() {
            @Override
            public void setValue(double newValue) {
                int idx = (int) Math.round(newValue);
                canvas.setConfigIndex(idx);
                colorSchemeDeviation.setConfigIndex(idx);
                colorSchemeDirection.setConfigIndex(idx);
                dsPrevTime.setPrevConfigIndex(idx);
                atomFilterDeviation.setConfigIndex(idx);
            }

            @Override
            public double getValue() {
                return colorSchemeDeviation.getConfigIndex();
            }

            @Override
            public Dimension getDimension() {
                return Null.DIMENSION;
            }

            @Override
            public String getLabel() {
                return "previous config";
            }
        });

        prevConfigSlider.setNMajor(5);
        prevConfigSlider.setMaximum(30);
        prevConfigSlider.setShowValues(true);
        prevConfigSlider.setEditValues(true);

        DeviceCheckBox linearCheckbox = new DeviceCheckBox("linear", new ModifierBoolean() {
            @Override
            public void setBoolean(boolean b) {
                if (b) {
                    canvas.setConfigStorage(configStorageLinear);
                    colorSchemeDeviation.setConfigStorage(configStorageLinear);
                    colorSchemeDirection.setConfigStorage(configStorageLinear);
                    dsPrevTime.setConfigStorage(configStorageLinear);
                    atomFilterDeviation.setConfigStorage(configStorageLinear);
                    prevConfigSlider.setMaximum(1000);
                }else{
                    canvas.setConfigStorage(configStorage);
                    colorSchemeDeviation.setConfigStorage(configStorage);
                    colorSchemeDirection.setConfigStorage(configStorage);
                    dsPrevTime.setConfigStorage(configStorage);
                    atomFilterDeviation.setConfigStorage(configStorage);
                    prevConfigSlider.setMaximum(30);
                }
                dbox.repaint();
            }

            @Override
            public boolean getBoolean() {
                return dsPrevTime.getConfigStorage() == configStorageLinear;
            }
        });
        linearCheckbox.setController(sim.getController());
        add(linearCheckbox);

        prevConfigSlider.setShowBorder(true);
        add(prevConfigSlider);



        DisplayTextBox displayPrevTime = new DisplayTextBox();
        DataPumpListener pumpPrevTime = new DataPumpListener(dsPrevTime, displayPrevTime, 1);
        sim.integrator.getEventManager().addListener(pumpPrevTime);
        add(displayPrevTime);

        DeviceSlider corIntervalSlider = new DeviceSlider(sim.getController(),null);
        corIntervalSlider.setShowBorder(true);
        corIntervalSlider.setNMajor(5);
        corIntervalSlider.setMaximum(20);
        corIntervalSlider.setMinimum(0);
        corIntervalSlider.setLabel("log2(sample interval (steps))");
        add(corIntervalSlider);


        DeviceSlider filterSlider = new DeviceSlider(sim.getController(), new Modifier() {
            @Override
            public void setValue(double newValue) {
                atomFilterDeviation.setMinDistance(newValue);
                dbox.repaint();
            }

            @Override
            public double getValue() {
                return atomFilterDeviation.getMinDistance();
            }

            @Override
            public Dimension getDimension() {
                return Length.DIMENSION;
            }

            @Override
            public String getLabel() {
                return "filter minDistance";
            }
        });
        filterSlider.setMaximum(2);
        filterSlider.setShowBorder(true);
        filterSlider.setPrecision(1);
        filterSlider.setNMajor(5);
        add(filterSlider);




        DeviceSlider stringSlider = new DeviceSlider(sim.getController(), new Modifier() {
            @Override
            public void setValue(double newValue) {
                colorSchemeDeviation.setDrString(newValue);
                dbox.repaint();
            }

            @Override
            public double getValue() {
                return colorSchemeDeviation.getDrString();
            }

            @Override
            public Dimension getDimension() {
                return Length.DIMENSION;
            }

            @Override
            public String getLabel() {
                return "filter strings";
            }
        });
        stringSlider.setPrecision(1);
        stringSlider.setMaximum(0.8);
        stringSlider.setShowBorder(true);
        stringSlider.setNMajor(4);
        add(stringSlider);


        DeviceCheckBox colorCheckboxStrings = new DeviceCheckBox("color strings", new ModifierBoolean() {
            @Override
            public void setBoolean(boolean b) {
                colorSchemeDeviation.setIsString(b);
                dbox.repaint();
            }

            @Override
            public boolean getBoolean() {
                return colorSchemeDeviation.getIsString();
            }
        });
        colorCheckboxStrings.setController(sim.getController());
        add(colorCheckboxStrings);



        DeviceCheckBox colorCheckbox = new DeviceCheckBox("color by displacement", null);
        colorCheckbox.setController(sim.getController());
        add(colorCheckbox);
        DeviceCheckBox colorDirectionCheckbox = new DeviceCheckBox("color by direction", new ModifierBoolean() {
            @Override
            public void setBoolean(boolean b) {
                if ((dbox.getColorScheme() == colorSchemeDirection) == b) return;
                if (b) {
                    if (sim.integrator.isIsothermal()) return;
                    colorCheckbox.setState(false);
                    dbox.setColorScheme(colorSchemeDirection);
                } else {
                    dbox.setColorScheme(colorScheme);
                }
                dbox.repaint();
            }

            @Override
            public boolean getBoolean() {
                return dbox.getColorScheme() == colorSchemeDirection;
            }
        });
        colorCheckbox.setModifier(new ModifierBoolean() {
            @Override
            public void setBoolean(boolean b) {
                if ((dbox.getColorScheme() == colorSchemeDeviation) == b) return;
                if (b) {
                    if (sim.integrator.isIsothermal()) return;
                    colorDirectionCheckbox.setState(false);
                    dbox.setColorScheme(colorSchemeDeviation);
                } else {
                    dbox.setColorScheme(colorScheme);
                }
                dbox.repaint();
            }

            @Override
            public boolean getBoolean() {
                return dbox.getColorScheme() == colorSchemeDeviation;
            }
        });
        colorCheckbox.setController(sim.getController());
        add(colorCheckbox);
        colorDirectionCheckbox.setController(sim.getController());
        add(colorDirectionCheckbox);
        if (sim.getSpace().getD() == 3) {
            JPanel xyzPanel = new JPanel(new GridBagLayout());
            GridBagConstraints h = SimulationPanel.getHorizGBC();
            String[] labels = new String[]{"x", "y", "z"};
            DeviceCheckBox[] axisCheckbox = new DeviceCheckBox[3];
            for (int i = 0; i < labels.length; i++) {
                axisCheckbox[i] = new DeviceCheckBox(labels[i], null);
            }
            for (int i = 0; i < labels.length; i++) {
                final int k = i;
                axisCheckbox[k].setModifier(new ModifierBoolean() {
                    @Override
                    public void setBoolean(boolean b) {
                        if (!b) {
                            if (colorSchemeDirection.getAxis() == k) throw new RuntimeException("oops");
                            return;
                        }
                        colorSchemeDirection.setAxis(k);
                        for (int j = 0; j < 3; j++) {
                            if (k == j) continue;
                            axisCheckbox[j].setState(false);
                        }
                        dbox.repaint();
                    }

                    @Override
                    public boolean getBoolean() {
                        return colorSchemeDirection.getAxis() == k;
                    }
                });
                xyzPanel.add(axisCheckbox[k].graphic());
            }
            getPanel().controlPanel.add(xyzPanel, SimulationPanel.getVertGBC());
        }

        DeviceCheckBox showDispCheckbox = new DeviceCheckBox("show displacement", new ModifierBoolean() {
            @Override
            public void setBoolean(boolean b) {
                if (canvas.getDrawDisplacement() == b || sim.integrator.isIsothermal()) return;
                canvas.setDrawDisplacement(b);
                dbox.repaint();
                diameterHash.setFac(b ? 0.5 : 1.0);
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
        repaintAction2.setInterval(REPAINT_INTERVAL);
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
        AccumulatorAutocorrelationPTensor dpAutocor = new AccumulatorAutocorrelationPTensor(256, sim.integrator.getTimeStep());
        pTensorFork.addDataSink(dpAutocor);
        dpAutocor.setPushInterval(16384);
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
        DataProcessorErrorBar pAutoCorErr = new DataProcessorErrorBar("err+");
        dpAutocor.getAvgErrFork().addDataSink(pAutoCorErr);


        //Gs: total
        int gsUpdateInterval = sim.getSpace().D() == 2 ? 10000 : 2000;
        double xGsMax = 2;
        MeterGs meterGs = new MeterGs(configStorageLinear);
        meterGs.setPrevSampleIndex(128);
        configStorageLinear.addListener(meterGs);
        meterGs.getXDataSource().setXMax(xGsMax);
        meterGs.reset();
        DisplayPlot gsPlot = new DisplayPlot();
        DataPumpListener pumpGs = new DataPumpListener(meterGs, gsPlot.getDataSet().makeDataSink(), gsUpdateInterval);
        sim.integrator.getEventManager().addListener(pumpGs);
        gsPlot.setLabel(" Gs ");
        add(gsPlot);
        //Gs: A
        MeterGs meterGsA = new MeterGs(configStorageLinear);
        meterGsA.setPrevSampleIndex(128);
        meterGsA.setAtomTypes(sim.speciesA.getLeafType());
        configStorageLinear.addListener(meterGsA);
        meterGsA.getXDataSource().setXMax(xGsMax);
        meterGsA.reset();
        DataPumpListener pumpGsA = new DataPumpListener(meterGsA, gsPlot.getDataSet().makeDataSink(), gsUpdateInterval);
        sim.integrator.getEventManager().addListener(pumpGsA);

        //Gs: B
        MeterGs meterGsB = new MeterGs(configStorageLinear);
        meterGsB.setPrevSampleIndex(128);
        meterGsB.setAtomTypes(sim.speciesB.getLeafType());
        configStorageLinear.addListener(meterGsB);
        meterGsB.getXDataSource().setXMax(xGsMax);
        meterGsB.reset();
        DataPumpListener pumpGsB = new DataPumpListener(meterGsB, gsPlot.getDataSet().makeDataSink(), gsUpdateInterval);
        sim.integrator.getEventManager().addListener(pumpGsB);

        gsPlot.setLegend(new DataTag[]{meterGs.getTag()}, "total");
        gsPlot.setLegend(new DataTag[]{meterGsA.getTag()}, "A");
        gsPlot.setLegend(new DataTag[]{meterGsB.getTag()}, "B");


        //Gs slider
        DeviceSlider gsPrevSampleSlider = new DeviceSlider(sim.getController(), new Modifier() {
            @Override
            public void setValue(double newValue) {
                if (newValue == getValue()) return;
                int log2prevConfig = (int) Math.round(newValue);
                int prevConfig = 1 << log2prevConfig;
                meterGs.setPrevSampleIndex(prevConfig);
                meterGsA.setPrevSampleIndex(prevConfig);
                meterGsB.setPrevSampleIndex(prevConfig);
            }

            @Override
            public double getValue() {
                int prevConfig = meterGs.getPrevSampleIndex();
                for (int i = 0; i <= 30; i++) {
                    if (1 << i >= prevConfig) return i;
                }
                throw new RuntimeException("oops");
            }

            @Override
            public Dimension getDimension() {
                return Null.DIMENSION;
            }

            @Override
            public String getLabel() {
                return null;
            }
        });
        gsPrevSampleSlider.setShowBorder(true);
        gsPrevSampleSlider.setShowValues(true);
        gsPrevSampleSlider.setNMajor(5);
        gsPrevSampleSlider.setMaximum(20);
        gsPrevSampleSlider.setMinimum(0);
        gsPrevSampleSlider.setLabel("log2(previous sample)");


        JPanel gsPanel = (JPanel) gsPlot.graphic();
        gsPanel.setLayout(new GridBagLayout());
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridheight = 1;
        gsPanel.add(gsPlot.getPlot(), gbc);
        gbc.insets = new Insets(20, 0, 0, 0);
        gbc.gridheight = 1;
        gbc.gridx = 1;
        gsPanel.add(gsPrevSampleSlider.graphic(), gbc);
        gbc.gridy = 1;
        gbc.gridx = gbc.gridy = 0;
        gbc.gridheight = 1;
        gbc.insets = new Insets(0, 0, 0, 0);



        int corUpdateInterval = sim.getSpace().D() == 2 ? 10000 : 2000;
        double xCorMax = 5;

        MeterCorrelation meterCorrelation = new MeterCorrelation(configStorageLinear);
        meterCorrelation.setPrevSampleIndex(128);
        configStorageLinear.addListener(meterCorrelation);
        meterCorrelation.getXDataSource().setXMax(xCorMax);
        meterCorrelation.reset();
        DisplayPlot correlationPlot = new DisplayPlot();
        DataPumpListener pumpCorrelation = new DataPumpListener(meterCorrelation, correlationPlot.getDataSet().makeDataSink(), corUpdateInterval);
        sim.integrator.getEventManager().addListener(pumpCorrelation);
        correlationPlot.setLabel("cor");
        add(correlationPlot);

        MeterCorrelation meterCorrelationAA = new MeterCorrelation(configStorageLinear);
        meterCorrelationAA.setAtomTypes(sim.speciesA.getLeafType(), sim.speciesA.getLeafType());
        meterCorrelationAA.setPrevSampleIndex(128);
        configStorageLinear.addListener(meterCorrelationAA);
        meterCorrelationAA.getXDataSource().setXMax(xCorMax);
        meterCorrelationAA.reset();
        DataPumpListener pumpCorrelationAA = new DataPumpListener(meterCorrelationAA, correlationPlot.getDataSet().makeDataSink(), corUpdateInterval);
        sim.integrator.getEventManager().addListener(pumpCorrelationAA);
        MeterCorrelation meterCorrelationAB = new MeterCorrelation(configStorageLinear);
        meterCorrelationAB.setAtomTypes(sim.speciesA.getLeafType(), sim.speciesB.getLeafType());
        meterCorrelationAB.setPrevSampleIndex(128);
        configStorageLinear.addListener(meterCorrelationAB);
        meterCorrelationAB.getXDataSource().setXMax(xCorMax);
        meterCorrelationAB.reset();
        DataPumpListener pumpCorrelationAB = new DataPumpListener(meterCorrelationAB, correlationPlot.getDataSet().makeDataSink(), corUpdateInterval);
        sim.integrator.getEventManager().addListener(pumpCorrelationAB);
        MeterCorrelation meterCorrelationBB = new MeterCorrelation(configStorageLinear);
        meterCorrelationBB.setAtomTypes(sim.speciesB.getLeafType(), sim.speciesB.getLeafType());
        meterCorrelationBB.setPrevSampleIndex(128);
        configStorageLinear.addListener(meterCorrelationBB);
        meterCorrelationBB.getXDataSource().setXMax(xCorMax);
        meterCorrelationBB.reset();
        DataPumpListener pumpCorrelationBB = new DataPumpListener(meterCorrelationBB, correlationPlot.getDataSet().makeDataSink(), corUpdateInterval);
        sim.integrator.getEventManager().addListener(pumpCorrelationBB);

        correlationPlot.setLegend(new DataTag[]{meterCorrelation.getTag()}, "total");
        correlationPlot.setLegend(new DataTag[]{meterCorrelationAA.getTag()}, "AA");
        correlationPlot.setLegend(new DataTag[]{meterCorrelationAB.getTag()}, "AB");
        correlationPlot.setLegend(new DataTag[]{meterCorrelationBB.getTag()}, "BB");

        MeterCorrelation meterCorrelationPerp = new MeterCorrelation(configStorageLinear, MeterCorrelation.CorrelationType.PERPENDICULAR);
        meterCorrelationPerp.setPrevSampleIndex(128);
        configStorageLinear.addListener(meterCorrelationPerp);
        meterCorrelationPerp.getXDataSource().setXMax(xCorMax);
        meterCorrelationPerp.reset();
        DataPumpListener pumpCorrelationPerp = new DataPumpListener(meterCorrelationPerp, correlationPlot.getDataSet().makeDataSink(), corUpdateInterval);
        sim.integrator.getEventManager().addListener(pumpCorrelationPerp);
        correlationPlot.setLegend(new DataTag[]{meterCorrelationPerp.getTag()}, "_|_");
        MeterCorrelation meterCorrelationPar = new MeterCorrelation(configStorageLinear, MeterCorrelation.CorrelationType.PARALLEL);
        meterCorrelationPar.setPrevSampleIndex(128);
        configStorageLinear.addListener(meterCorrelationPar);
        meterCorrelationPar.getXDataSource().setXMax(xCorMax);
        meterCorrelationPar.reset();
        DataPumpListener pumpCorrelationPar = new DataPumpListener(meterCorrelationPar, correlationPlot.getDataSet().makeDataSink(), corUpdateInterval);
        sim.integrator.getEventManager().addListener(pumpCorrelationPar);
        correlationPlot.setLegend(new DataTag[]{meterCorrelationPar.getTag()}, "||");

        MeterCorrelation meterCorrelationMag = new MeterCorrelation(configStorageLinear, MeterCorrelation.CorrelationType.MAGNITUDE);
        meterCorrelationMag.setPrevSampleIndex(128);
        configStorageLinear.addListener(meterCorrelationMag);
        meterCorrelationMag.getXDataSource().setXMax(xCorMax);
        meterCorrelationMag.reset();
        DataPumpListener pumpCorrelationMag = new DataPumpListener(meterCorrelationMag, correlationPlot.getDataSet().makeDataSink(), corUpdateInterval);
        sim.integrator.getEventManager().addListener(pumpCorrelationMag);
        correlationPlot.setLegend(new DataTag[]{meterCorrelationMag.getTag()}, "|r|");


        DeviceSlider corPrevSampleSlider = new DeviceSlider(sim.getController(), new Modifier() {
            @Override
            public void setValue(double newValue) {
                if (newValue == getValue()) return;
                int log2prevConfig = (int) Math.round(newValue);
                int prevConfig = 1 << log2prevConfig;
                meterCorrelationAA.setPrevSampleIndex(prevConfig);
                meterCorrelationAB.setPrevSampleIndex(prevConfig);
                meterCorrelationBB.setPrevSampleIndex(prevConfig);
                meterCorrelation.setPrevSampleIndex(prevConfig);
                meterCorrelationPerp.setPrevSampleIndex(prevConfig);
                meterCorrelationPar.setPrevSampleIndex(prevConfig);
                meterCorrelationMag.setPrevSampleIndex(prevConfig);
            }

            @Override
            public double getValue() {
                int prevConfig = meterCorrelation.getPrevSampleIndex();
                for (int i = 0; i <= 30; i++) {
                    if (1 << i >= prevConfig) return i;
                }
                throw new RuntimeException("oops");
            }

            @Override
            public Dimension getDimension() {
                return Null.DIMENSION;
            }

            @Override
            public String getLabel() {
                return null;
            }
        });
        corPrevSampleSlider.setShowBorder(true);
        corPrevSampleSlider.setShowValues(true);
        corPrevSampleSlider.setNMajor(5);
        corPrevSampleSlider.setMaximum(10);
        corPrevSampleSlider.setMinimum(0);
        corPrevSampleSlider.setLabel("log2(previous sample)");

        corIntervalSlider.setModifier(new Modifier() {
            @Override
            public void setValue(double newValue) {
                int log2interval = (int) Math.round(newValue);
                int interval = 1 << log2interval;
                configStorageLinear.setSampleInterval(interval);
                meterCorrelationAA.reset();
                meterCorrelationAB.reset();
                meterCorrelationBB.reset();
                meterCorrelation.reset();
                meterCorrelationPerp.reset();
                meterCorrelationPar.reset();
                meterCorrelationMag.reset();
                meterGs.reset();
                meterGsA.reset();
                meterGsB.reset();
            }

            @Override
            public double getValue() {
                int interval = configStorageLinear.getSampleInterval();
                for (int i = 0; i <= 30; i++) {
                    if (1 << i >= interval) return i;
                }
                throw new RuntimeException("oops");
            }

            @Override
            public Dimension getDimension() {
                return Null.DIMENSION;
            }

            @Override
            public String getLabel() {
                return "log2(sample interval (steps))";
            }
        });

        JPanel corPanel = (JPanel) correlationPlot.graphic();
        corPanel.setLayout(new GridBagLayout());
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridheight = 2;
        corPanel.add(correlationPlot.getPlot(), gbc);
        gbc.insets = new Insets(20, 0, 0, 0);
        gbc.gridheight = 1;
        gbc.gridx = 1;
        corPanel.add(corPrevSampleSlider.graphic(), gbc);
        gbc.gridx = gbc.gridy = 0;
        gbc.gridheight = 1;
        gbc.insets = new Insets(0, 0, 0, 0);


        DisplayPlot plotPTensorAutocor = new DisplayPlot();
        plotPTensorAutocor.setLabel("P Tensor autocor");
        plotPTensorAutocor.setLegend(new DataTag[]{dpAutocor.getTag()}, "avg");
        dpAutocor.addDataSink(plotPTensorAutocor.getDataSet().makeDataSink());
        add(plotPTensorAutocor);
        DeviceSlider nMaxSlider = new DeviceSlider(sim.getController(), new Modifier() {
            @Override
            public void setValue(double newValue) {
                int log2nMax = (int) Math.round(newValue);
                int nMax = 1 << log2nMax;
                dpAutocor.setNMax(nMax);
            }

            @Override
            public double getValue() {
                int nMax = dpAutocor.getNMax();
                for (int i = 0; i <= 30; i++) {
                    if (1 << i >= nMax) return i;
                }
                throw new RuntimeException("oops");
            }

            @Override
            public Dimension getDimension() {
                return Null.DIMENSION;
            }

            @Override
            public String getLabel() {
                return null;
            }
        });
        nMaxSlider.setShowBorder(true);
        nMaxSlider.setShowValues(true);
        nMaxSlider.setNMajor(5);
        nMaxSlider.setMaximum(30);
        nMaxSlider.setMinimum(0);
        nMaxSlider.setLabel("log2(Max lag (steps))");
        DeviceSlider pushIntervalSlider = new DeviceSlider(sim.getController(), new Modifier() {
            @Override
            public void setValue(double newValue) {
                int log2pi = (int) Math.round(newValue);
                int pi = 1 << log2pi;
                dpAutocor.setPushInterval(pi);
            }

            @Override
            public double getValue() {
                long pi = dpAutocor.getPushInterval();
                for (int i = 0; i <= 30; i++) {
                    if (1 << i >= pi) return i;
                }
                throw new RuntimeException("oops");
            }

            @Override
            public Dimension getDimension() {
                return Null.DIMENSION;
            }

            @Override
            public String getLabel() {
                return null;
            }
        });
        pushIntervalSlider.setShowBorder(true);
        pushIntervalSlider.setShowValues(true);
        pushIntervalSlider.setNMajor(5);
        pushIntervalSlider.setMaximum(30);
        pushIntervalSlider.setMinimum(0);
        pushIntervalSlider.setLabel("log2(Push interval (steps))");
        JPanel ptacPanel = (JPanel) plotPTensorAutocor.graphic();
        ptacPanel.remove(plotPTensorAutocor.getPlot());
        ptacPanel.setLayout(new GridBagLayout());
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridheight = 2;
        ptacPanel.add(plotPTensorAutocor.getPlot(), gbc);
        gbc.insets = new Insets(20, 0, 0, 0);
        gbc.gridheight = 1;
        gbc.gridx = 1;
        ptacPanel.add(nMaxSlider.graphic(), gbc);
        gbc.gridy = 1;
        ptacPanel.add(pushIntervalSlider.graphic(), gbc);

        pAutoCorErr.setDataSink(plotPTensorAutocor.getDataSet().makeDataSink());
        plotPTensorAutocor.setLegend(new DataTag[]{dpAutocor.getAvgErrFork().getTag(), pAutoCorErr.getTag()}, "err+");

        DisplayTextBoxesCAE peDisplay = null;
        if (sim.potentialChoice != SimGlass.PotentialChoice.HS) {
            AccumulatorHistory historyPE = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyPE.setTimeDataSource(timeCounter);
            peFork.addDataSink(historyPE);
            DisplayPlot plotPE = new DisplayPlot();
            historyPE.addDataSink(plotPE.getDataSet().makeDataSink());
            plotPE.setLabel("PE");
            plotPE.setUnit(peUnit);
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
        plotMSD.setLegend(new DataTag[]{meterMSD.getTag()}, "total");
        sim.integrator.getEventManager().addListener(pumpMSD);
        plotMSD.setLabel("MSD");
        plotMSD.getPlot().setYLog(true);
        plotMSD.getPlot().setXLog(true);
        add(plotMSD);
        DataSourceMSD meterMSDA = new DataSourceMSD(configStorageMSD, sim.speciesA.getLeafType());
        configStorageMSD.addListener(meterMSDA);
        DataPumpListener pumpMSDA = new DataPumpListener(meterMSDA, plotMSD.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpMSDA);
        plotMSD.setLegend(new DataTag[]{meterMSDA.getTag()}, "A");
        DataSourceMSD meterMSDB = new DataSourceMSD(configStorageMSD, sim.speciesB.getLeafType());
        configStorageMSD.addListener(meterMSDB);
        DataPumpListener pumpMSDB = new DataPumpListener(meterMSDB, plotMSD.getDataSet().makeDataSink(), 1000);
        plotMSD.setLegend(new DataTag[]{meterMSDB.getTag()}, "B");
        sim.integrator.getEventManager().addListener(pumpMSDB);


        DataSourceAlpha2 meterAlpha2 = new DataSourceAlpha2(configStorageMSD);
        configStorageMSD.addListener(meterAlpha2);
        DisplayPlot plotAlpha2 = new DisplayPlot();
        DataPumpListener pumpAlpha2 = new DataPumpListener(meterAlpha2, plotAlpha2.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpAlpha2);
        plotAlpha2.setLabel("alpha2");
        plotAlpha2.getPlot().setXLog(true);
        add(plotAlpha2);


        //Fs: TOTAL
        DataSourceFs meterFs = new DataSourceFs(configStorageMSD);
        configStorageMSD.addListener(meterFs);
        DisplayPlot plotFs = new DisplayPlot();
        DataPumpListener pumpFs = new DataPumpListener(meterFs, plotFs.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpFs);
        plotFs.setLabel("Fs");
        plotFs.getPlot().setXLog(true);
        add(plotFs);

        //Fs: A
        DataSourceFs meterFsA = new DataSourceFs(configStorageMSD);
        meterFsA.setAtomType(sim.speciesA.getLeafType());
        configStorageMSD.addListener(meterFsA);
        DataPumpListener pumpFsA = new DataPumpListener(meterFsA, plotFs.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpFsA);

        //Fs: B
        DataSourceFs meterFsB = new DataSourceFs(configStorageMSD);
        meterFsB.setAtomType(sim.speciesB.getLeafType());
        configStorageMSD.addListener(meterFsB);
        DataPumpListener pumpFsB = new DataPumpListener(meterFsB, plotFs.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpFsB);


        plotFs.setLegend(new DataTag[]{meterFs.getTag()}, "total");
        plotFs.setLegend(new DataTag[]{meterFsA.getTag()}, "A");
        plotFs.setLegend(new DataTag[]{meterFsB.getTag()}, "B");


        //************* Lay out components ****************//

        DeviceCheckBox swapCheckbox = new DeviceCheckBox("isothermal", new ModifierBoolean() {
            @Override
            public void setBoolean(boolean b) {

                if (sim.integrator.isIsothermal() == b) return;
                if (b) {
                    dbox.setAtomTestDoDisplay(null);
                    sim.integrator.setIsothermal(true);
                    sim.integrator.setIntegratorMC(sim.integratorMC, 10000);
                    dbox.setColorScheme(colorScheme);
                    canvas.setDrawDisplacement(false);
                    configStorageLinear.reset();
                    configStorageLinear.setEnabled(false);
                    configStorage.reset();
                    configStorage.setEnabled(false);
                    configStorageMSD.reset();
                    configStorageMSD.setEnabled(false);
                    dpAutocor.reset();
                    meterCorrelationAA.reset();
                    meterCorrelationAB.reset();
                    meterCorrelationBB.reset();
                    meterCorrelation.reset();
                    meterCorrelationPerp.reset();
                    meterCorrelationPar.reset();
                    meterCorrelationMag.reset();
                    diameterHash.setFac(1.0);
                } else {
                    dbox.setAtomTestDoDisplay(atomFilterDeviation);
                    sim.integrator.setIntegratorMC(null, 0);
                    sim.integrator.setIsothermal(false);
                    configStorageLinear.reset();
                    configStorageLinear.setEnabled(true);
                    configStorage.reset();
                    configStorage.setEnabled(true);
                    configStorageMSD.reset();
                    configStorageMSD.setEnabled(true);
                    if (colorCheckbox.getState()) dbox.setColorScheme(colorSchemeDeviation);
                    else if (colorDirectionCheckbox.getState()) dbox.setColorScheme(colorSchemeDirection);
                    if (showDispCheckbox.getState()) canvas.setDrawDisplacement(true);
                    dpAutocor.reset();
                    meterCorrelationAA.reset();
                    meterCorrelationAB.reset();
                    meterCorrelationBB.reset();
                    meterCorrelation.reset();
                    meterCorrelationPerp.reset();
                    meterCorrelationPar.reset();
                    meterCorrelationMag.reset();
                    diameterHash.setFac(showDispCheckbox.getState() ? 0.5 : 1.0);
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
        configStorageLinear.reset();
        configStorageMSD.reset();
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
            params.potential = SimGlass.PotentialChoice.LJ;
            params.nA = 800;
            params.nB = 200;
            params.density = 1.25; // 3D 1.25;
            params.D = 3;
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