/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.glass;

import etomica.action.IAction;
import etomica.atom.AtomType;
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
import etomica.modifier.ModifierGeneral;
import etomica.space.Vector;
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

    private final static String APP_NAME = " Molecular Dynamics";
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

        ConfigurationStorage configStorageMSDPerc = new ConfigurationStorage(sim.box, ConfigurationStorage.StorageType.MSD);
        configStorageMSDPerc.setEnabled(false); // start isothermal
        sim.integrator.getEventManager().addListener(configStorageMSDPerc);

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
        diameterHash.setDiameter(sim.speciesB.getLeafType(), sim.sigmaB);
        diameterHash.setDiameter(sim.speciesA.getLeafType(), 1);
        dbox.setDiameterHash(diameterHash);

        AtomTestDeviation atomFilterDeviation = new AtomTestDeviation(sim.box, configStorage);
        AtomTestDeviation atomFilterDeviationPerc = new AtomTestDeviation(sim.box, configStorageMSDPerc);

        ColorSchemeDeviation colorSchemeDeviation = new ColorSchemeDeviation(sim.box, configStorage);
        ColorSchemeDirection colorSchemeDirection = new ColorSchemeDirection(sim.box, configStorage);
        ColorSchemeCluster colorSchemeCluster = new ColorSchemeCluster(sim.box, atomFilterDeviation);

        DataSourcePrevTime dsPrevTime = new DataSourcePrevTime(configStorage);
        DisplayTextBox displayPrevTime = new DisplayTextBox();
        DataPumpListener pumpPrevTime = new DataPumpListener(dsPrevTime, displayPrevTime, 1);
        sim.integrator.getEventManager().addListener(pumpPrevTime);

        DeviceSlider prevConfigSlider = new DeviceSlider(sim.getController(), new Modifier() {
            @Override
            public void setValue(double newValue) {
                int idx = (int) Math.round(newValue);
                canvas.setConfigIndex(idx);
                colorSchemeDeviation.setConfigIndex(idx);
                colorSchemeDirection.setConfigIndex(idx);
                dsPrevTime.setPrevConfigIndex(idx);
                atomFilterDeviation.setConfigIndex(idx);
                pumpPrevTime.actionPerformed();
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
                pumpPrevTime.actionPerformed();
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
        filterSlider.setValue(0.4);
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

        DeviceCheckBox immobileCheckbox = new DeviceCheckBox("show immobile", null);
        immobileCheckbox.setController(sim.getController());
        add(immobileCheckbox);

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
                    immobileCheckbox.setState(false);
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
                    immobileCheckbox.setState(false);
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
        immobileCheckbox.setModifier(new ModifierBoolean() {
            @Override
            public void setBoolean(boolean b) {
                if ((dbox.getColorScheme() == colorSchemeCluster) == b) return;
                atomFilterDeviation.setDoMobileOnly(!b);
                if (b) {
                    if (sim.integrator.isIsothermal()) return;
                    colorDirectionCheckbox.setState(false);
                    colorCheckbox.setState(false);
                    dbox.setColorScheme(colorSchemeCluster);
                } else {
                    dbox.setColorScheme(colorScheme);
                }
                dbox.repaint();
            }

            @Override
            public boolean getBoolean() {
                return !atomFilterDeviation.getDoMobileOnly();
            }
        });

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

        DeviceCheckBox flipDispCheckbox = new DeviceCheckBox("flip displacement", new ModifierBoolean() {
            @Override
            public void setBoolean(boolean b) {
                if (canvas.getFlipDisplacement() == b || sim.integrator.isIsothermal()) return;
                canvas.setFlipDisplacement(b);
                dbox.repaint();
            }

            @Override
            public boolean getBoolean() {
                return canvas.getFlipDisplacement();
            }
        });

        flipDispCheckbox.setController(sim.getController());
        add(flipDispCheckbox);

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

        AccumulatorPTensor pTensorAccum = new AccumulatorPTensor(sim.integrator, sim.integrator.getTimeStep());

        AccumulatorHistory energyHistory, peHistory;
        final AccumulatorAverageCollapsing peAccumulator;
        DisplayPlot ePlot, tPlot;
        DataFork peFork = null;
        DataFork tFork = null;
        int log2peInterval = 6;
        if (sim.potentialChoice != SimGlass.PotentialChoice.HS) {
            MeterEnergy eMeter = new MeterEnergy(sim.integrator.getPotentialMaster(), sim.box);
            energyHistory = new AccumulatorHistory();
            energyHistory.setTimeDataSource(timeCounter);
            DataPumpListener energyPump = new DataPumpListener(eMeter, energyHistory, 1 << log2peInterval);
            sim.integrator.getEventManager().addListener(energyPump);
            energyHistory.setPushInterval(5);
            dataStreamPumps.add(energyPump);

            MeterPotentialEnergy peMeter = new MeterPotentialEnergy(sim.integrator.getPotentialMaster(), sim.box);
            peHistory = new AccumulatorHistory();
            peHistory.setTimeDataSource(timeCounter);
            peAccumulator = new AccumulatorAverageCollapsing();
            peAccumulator.setPushInterval(10);
            peFork = new DataFork(new IDataSink[]{peHistory, peAccumulator});
            DataPumpListener pePump = new DataPumpListener(peMeter, peFork, 1 << log2peInterval);
            sim.integrator.getEventManager().addListener(pePump);
            peHistory.setPushInterval(5);
            dataStreamPumps.add(pePump);

            MeterKineticEnergy keMeter = new MeterKineticEnergy(sim.box);
            AccumulatorHistory keHistory = new AccumulatorHistory();
            keHistory.setTimeDataSource(timeCounter);
            DataFork keFork = new DataFork();
            DataPumpListener kePump = new DataPumpListener(keMeter, keFork, 1 << log2peInterval);
            keFork.addDataSink(keHistory);
            final AccumulatorAverage keAvg = new AccumulatorAverageCollapsing();
            keFork.addDataSink(keAvg);
            sim.integrator.getEventManager().addListener(kePump);
            keHistory.setPushInterval(5);
            dataStreamPumps.add(kePump);

            ePlot = new DisplayPlot();
            energyHistory.setDataSink(ePlot.getDataSet().makeDataSink());
            ePlot.setLegend(new DataTag[]{energyHistory.getTag()}, "Total");
            peHistory.setDataSink(ePlot.getDataSet().makeDataSink());
            ePlot.setLegend(new DataTag[]{peHistory.getTag()}, "Potential");
            keHistory.setDataSink(ePlot.getDataSet().makeDataSink());
            ePlot.setLegend(new DataTag[]{keHistory.getTag()}, "Kinetic");

            ePlot.getPlot().setTitle("Energy History");
            ePlot.setDoLegend(true);
            ePlot.setLabel("Energy");


            //Temperature:
            MeterTemperature tMeter = new MeterTemperature(sim.box, sim.getSpace().getD());
            AccumulatorHistory tHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
            tHistory.setTimeDataSource(timeCounter);
            AccumulatorAverageCollapsing tAccumulator = new AccumulatorAverageCollapsing();
            tAccumulator.setPushInterval(10);
            tFork = new DataFork(new IDataSink[]{tHistory, tAccumulator});
            DataPumpListener tPump = new DataPumpListener(tMeter, tFork, 1 << log2peInterval);
            sim.integrator.getEventManager().addListener(tPump);
            tHistory.setPushInterval(5);
            dataStreamPumps.add(tPump);

            tPlot = new DisplayPlot();
            tHistory.setDataSink(tPlot.getDataSet().makeDataSink());
            tPlot.setDoLegend(false);
            tPlot.setLabel("T");

            tAccumulator.addDataSink(pTensorAccum.makeTemperatureSink(), new AccumulatorAverage.StatType[]{tAccumulator.AVERAGE});

        } else {
            peAccumulator = null;
            ePlot = null;
            tPlot = null;
        }

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

        //unnormalized AC of all stress tensor components
        AccumulatorAutocorrelationPTensor dpAutocor = new AccumulatorAutocorrelationPTensor(256, sim.integrator.getTimeStep());
        if (sim.box.getLeafList().size() > 200 && sim.potentialChoice != SimGlass.PotentialChoice.HS) {
            pTensorFork.addDataSink(dpAutocor);
        }

        //normalized AC of shear stress components
        AccumulatorAutocorrelationShearStress dpxyAutocor = new AccumulatorAutocorrelationShearStress(256, sim.integrator.getTimeStep());
        if (sim.box.getLeafList().size() > 200 && sim.potentialChoice != SimGlass.PotentialChoice.HS) {
            pTensorFork.addDataSink(dpxyAutocor);
        }


        DisplayPlot plotPTensorAccum = new DisplayPlot();
        plotPTensorAccum.setLabel("viscosity(t)");
        plotPTensorAccum.setDoLegend(false);
        plotPTensorAccum.getPlot().setXLog(true);
        plotPTensorAccum.getPlot().setYLog(true);
        add(plotPTensorAccum);
        DataPumpListener pTensorAccumPump = new DataPumpListener(pTensorAccum,plotPTensorAccum.getDataSet().makeDataSink(),1000);
        pTensorFork.addDataSink(pTensorAccum);
        sim.integrator.getEventManager().addListener(pTensorAccumPump);

        dpAutocor.setPushInterval(16384);
        dpxyAutocor.setPushInterval(16384);
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
        plotP.setDoLegend(false);
        add(plotP);
        DataProcessorErrorBar pAutoCorErr = new DataProcessorErrorBar("err+");
        dpAutocor.getAvgErrFork().addDataSink(pAutoCorErr);

        DataSourceMSDcorP dsMSDcorP = new DataSourceMSDcorP(sim.integrator);
        pFork.addDataSink(dsMSDcorP);

        DataSourceHisogram dsHistogramP = new DataSourceHisogram(sim.integrator);
        pFork.addDataSink(dsHistogramP);

        DataSourcePMSDHistory dsPMSDhistory = new DataSourcePMSDHistory(sim.integrator);
        pFork.addDataSink(dsPMSDhistory);

        DataSourceMSDcorP dsMSDcorU;
        if (sim.potentialChoice != SimGlass.PotentialChoice.HS && false) {
            // disabled.  U is not differently correlated than P
            dsMSDcorU = new DataSourceMSDcorP(sim.integrator, log2peInterval);
            peFork.addDataSink(dsMSDcorU);
        } else {
            dsMSDcorU = null;
        }



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


        //Strings
        DataSourceStrings meterStrings = new DataSourceStrings(configStorageMSD, 3, 30);
        configStorageMSD.addListener(meterStrings);
        DisplayPlot plotStrings = new DisplayPlot();
        DataPumpListener pumpStrings = new DataPumpListener(meterStrings, plotStrings.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpStrings);
        plotStrings.setLabel("strings");
        plotStrings.getPlot().setXLog(true);
        plotStrings.setDoLegend(false);
        add(plotStrings);

        //Percolation
        atomFilterDeviationPerc.setDoMobileOnly(false);
        DataSourcePercolation meterPerc = new DataSourcePercolation(configStorageMSDPerc, atomFilterDeviationPerc, 5, 30);
        configStorageMSDPerc.addListener(meterPerc);
        DisplayPlot plotPerc = new DisplayPlot();
        DataPumpListener pumpPerc = new DataPumpListener(meterPerc, plotPerc.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpPerc);

        DataPumpListener pumpImmFrac = new DataPumpListener(meterPerc.makeImmFractionSource(), plotPerc.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpImmFrac);


        plotPerc.setLegend(new DataTag[]{meterPerc.getTag()},"perc. prob.");
        plotPerc.setLegend(new DataTag[]{meterPerc.makeImmFractionSource().getTag()},"imm. frac.");
        plotPerc.setLabel("perc");
        plotPerc.getPlot().setXLog(true);
        add(plotPerc);

        //Percolation slider
        DeviceSlider percDrSlider = new DeviceSlider(sim.getController(), new Modifier() {
            @Override
            public void setValue(double newValue) {
                if (newValue == getValue()) return;
                configStorageMSDPerc.reset();
                atomFilterDeviationPerc.setMinDistance(newValue);
                meterPerc.reset(); // Needed ???
            }

            @Override
            public double getValue() {
                return atomFilterDeviationPerc.getMinDistance();
            }

            @Override
            public Dimension getDimension() {
                return Null.DIMENSION;
            }

            @Override
            public String getLabel() {
                return "filter perc. minDistance";
            }
        });
        percDrSlider.setShowBorder(true);
        percDrSlider.setShowValues(true);
        percDrSlider.setPrecision(1);
        percDrSlider.setNMajor(5);
        percDrSlider.setMaximum(2);
        percDrSlider.setMinimum(0);
        percDrSlider.setValue(0.4);
        percDrSlider.setLabel("immobile minDistance");

        //MinTime slider
        DeviceSlider percMinLog2StepSlider = new DeviceSlider(sim.getController(), new Modifier() {
            @Override
            public void setValue(double newValue) {
                if (newValue == getValue()) return;
                meterPerc.setLog2StepStart((int) newValue);
                configStorageMSDPerc.reset();
            }

            @Override
            public double getValue() {
                return meterPerc.getLog2StepStart();
            }

            @Override
            public Dimension getDimension() {
                return Null.DIMENSION;
            }

            @Override
            public String getLabel() {
                return "log2(min time (steps))";
            }
        });
        percMinLog2StepSlider.setShowBorder(true);
        percMinLog2StepSlider.setShowValues(true);
        percMinLog2StepSlider.setNMajor(6);
        percMinLog2StepSlider.setMaximum(30);
        percMinLog2StepSlider.setMinimum(0);
        percMinLog2StepSlider.setValue(5);
        percMinLog2StepSlider.setLabel("log2(min time (steps))");


        //MaxTime slider
        DeviceSlider percMaxLog2StepSlider = new DeviceSlider(sim.getController(), new Modifier() {
            @Override
            public void setValue(double newValue) {
                if (newValue == getValue()) return;
                meterPerc.setLog2StepEnd((int) newValue);
                configStorageMSDPerc.reset();
            }

            @Override
            public double getValue() {
                return meterPerc.getLog2StepEnd();
            }

            @Override
            public Dimension getDimension() {
                return Null.DIMENSION;
            }

            @Override
            public String getLabel() {
                return "log2(max time (steps))";
            }
        });
        percMaxLog2StepSlider.setShowBorder(true);
        percMaxLog2StepSlider.setShowValues(true);
        percMaxLog2StepSlider.setNMajor(6);
        percMaxLog2StepSlider.setMaximum(30);
        percMaxLog2StepSlider.setMinimum(0);
        percMaxLog2StepSlider.setValue(30);
        percMaxLog2StepSlider.setLabel("log2(max time (steps))");


        JPanel percPanel = (JPanel) plotPerc.graphic();
        percPanel.remove(plotPerc.getPlot());
        percPanel.setLayout(new GridBagLayout());
        GridBagConstraints gbcPerc = new GridBagConstraints();
        gbcPerc.gridx = 0;
        gbcPerc.gridy = 0;
        gbcPerc.gridheight = 3;
        percPanel.add(plotPerc.getPlot(), gbcPerc);
        gbcPerc.insets = new Insets(20, 0, 0, 0);
        gbcPerc.gridheight = 1;
        gbcPerc.gridx = 1;
        percPanel.add(percDrSlider.graphic(), gbcPerc);
        gbcPerc.gridy = 1;
        percPanel.add(percMinLog2StepSlider.getPanel(), gbcPerc);
        gbcPerc.gridy = 2;
        percPanel.add(percMaxLog2StepSlider.graphic(), gbcPerc);


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
                pumpPrevTime.actionPerformed();
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


        //unnormalized AC of all stress tensor components
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
        pushIntervalSlider.setValue(20);
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






        //normalized AC of shear stress components
        DisplayPlot plotPxyTensorAutocor = new DisplayPlot();
        plotPxyTensorAutocor.setLabel("Pxy autocor");
        plotPxyTensorAutocor.setLegend(new DataTag[]{dpxyAutocor.getTag()}, "avg");
        dpxyAutocor.addDataSink(plotPxyTensorAutocor.getDataSet().makeDataSink());
        add(plotPxyTensorAutocor);
        DeviceSlider nMaxSliderShear = new DeviceSlider(sim.getController(), new Modifier() {
            @Override
            public void setValue(double newValue) {
                int log2nMax = (int) Math.round(newValue);
                int nMax = 1 << log2nMax;
                dpxyAutocor.setNMax(nMax);
            }

            @Override
            public double getValue() {
                int nMax = dpxyAutocor.getNMax();
                for (int i = 0; i <= 30; i++) {
                    if (1 << i >= nMax) return i;
                }
                throw new RuntimeException("oops!");
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
        nMaxSliderShear.setShowBorder(true);
        nMaxSliderShear.setShowValues(true);
        nMaxSliderShear.setNMajor(5);
        nMaxSliderShear.setMaximum(30);
        nMaxSliderShear.setMinimum(0);
        nMaxSliderShear.setLabel("log2(Max lag (steps))");
        DeviceSlider pushIntervalSliderShear = new DeviceSlider(sim.getController(), new Modifier() {
            @Override
            public void setValue(double newValue) {
                int log2pi = (int) Math.round(newValue);
                int pi = 1 << log2pi;
                dpxyAutocor.setPushInterval(pi);
            }

            @Override
            public double getValue() {
                long pi = dpxyAutocor.getPushInterval();
                for (int i = 0; i <= 30; i++) {
                    if (1 << i >= pi) return i;
                }
                throw new RuntimeException("oops!");
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
        pushIntervalSliderShear.setShowBorder(true);
        pushIntervalSliderShear.setShowValues(true);
        pushIntervalSliderShear.setNMajor(5);
        pushIntervalSliderShear.setMaximum(30);
        pushIntervalSliderShear.setMinimum(0);
        pushIntervalSliderShear.setValue(10);
        pushIntervalSliderShear.setLabel("log2(Push interval (steps))");
        JPanel shearacPanel = (JPanel) plotPxyTensorAutocor.graphic();
        shearacPanel.remove(plotPxyTensorAutocor.getPlot());
        shearacPanel.setLayout(new GridBagLayout());
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridheight = 2;
        shearacPanel.add(plotPxyTensorAutocor.getPlot(), gbc);
        gbc.insets = new Insets(20, 0, 0, 0);
        gbc.gridheight = 1;
        gbc.gridx = 1;
        shearacPanel.add(nMaxSliderShear.graphic(), gbc);
        gbc.gridy = 1;
        shearacPanel.add(pushIntervalSliderShear.graphic(), gbc);







        //Potential energy
        DisplayTextBoxesCAE peDisplay = null;
        if (sim.potentialChoice != SimGlass.PotentialChoice.HS) {
            AccumulatorHistory historyPE = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyPE.setTimeDataSource(timeCounter);
            peFork.addDataSink(historyPE);
            DisplayPlot plotPE = new DisplayPlot();
            historyPE.addDataSink(plotPE.getDataSet().makeDataSink());
            plotPE.setLabel("PE");
            plotPE.setUnit(peUnit);
            plotPE.setDoLegend(false);
            add(plotPE);
            peDisplay = new DisplayTextBoxesCAE();
            peDisplay.setAccumulator(peAccumulator);
        }

        final DisplayTextBoxesCAE pDisplay = new DisplayTextBoxesCAE();
        pDisplay.setAccumulator(pAccumulator);




        //MSD
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
        plotAlpha2.setDoLegend(false);
        add(plotAlpha2);

        //FOld
        DataSourceFOld meterFOld = new DataSourceFOld(configStorageMSD);
        configStorageMSD.addListener(meterFOld);
        DisplayPlot plotF = new DisplayPlot();
        DataPumpListener pumpFOld = new DataPumpListener(meterFOld, plotF.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpFOld);


        //F - new
        DataSourceF meterF = new DataSourceF(configStorageMSD);
        configStorageMSD.addListener(meterF);
        DataPumpListener pumpF = new DataPumpListener(meterF, plotF.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpF);
        plotF.setLabel("F");
        plotF.getPlot().setXLog(true);
        add(plotF);


        plotF.setLegend(new DataTag[]{meterFOld.getTag()}, "old");
        plotF.setLegend(new DataTag[]{meterF.getTag()}, "new");



        //Fs: TOTAL
        double qx = 7.0;
        Vector qVec = sim.getSpace().makeVector();
        qVec.setX(0, qx);

        DataSourceFs meterFs = new DataSourceFs(configStorageMSD);
        meterFs.setQ(qVec);
        configStorageMSD.addListener(meterFs);
        DisplayPlot plotFs = new DisplayPlot();
        DataPumpListener pumpFs = new DataPumpListener(meterFs, plotFs.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpFs);
        plotFs.setLabel("Fs");
        plotFs.getPlot().setXLog(true);
        add(plotFs);

        //Fs: A
        DataSourceFs meterFsA = new DataSourceFs(configStorageMSD);
        meterFsA.setQ(qVec);
        meterFsA.setAtomType(sim.speciesA.getLeafType());
        configStorageMSD.addListener(meterFsA);
        DataPumpListener pumpFsA = new DataPumpListener(meterFsA, plotFs.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpFsA);

        //Fs: B
        DataSourceFs meterFsB = new DataSourceFs(configStorageMSD);
        meterFsB.setQ(qVec);
        meterFsB.setAtomType(sim.speciesB.getLeafType());
        configStorageMSD.addListener(meterFsB);
        DataPumpListener pumpFsB = new DataPumpListener(meterFsB, plotFs.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpFsB);


        plotFs.setLegend(new DataTag[]{meterFs.getTag()}, "total");
        plotFs.setLegend(new DataTag[]{meterFsA.getTag()}, "A");
        plotFs.setLegend(new DataTag[]{meterFsB.getTag()}, "B");




        DeviceSlider qSlider = new DeviceSlider(sim.getController(), new Modifier() {
            @Override
            public void setValue(double newValue) {
                Vector qVec = sim.getSpace().makeVector();
                qVec.setX(0, newValue);
                meterFs.setQ(qVec);
                meterFsA.setQ(qVec);
                meterFsB.setQ(qVec);
                configStorageMSD.reset();
                meterFs.reset();
                meterFsA.reset();
                meterFsB.reset();
                dbox.repaint();
            }

            @Override
            public double getValue() {
                return meterFs.getQ().getX(0);
            }

            @Override
            public Dimension getDimension() {
                return Length.DIMENSION;
            }

            @Override
            public String getLabel() {
                return "qx";
            }
        });
        qSlider.setMaximum(20);
        qSlider.setShowBorder(true);
        qSlider.setShowValues(true);
        qSlider.setEditValues(true);
        qSlider.setPrecision(1);
        qSlider.setNMajor(5);
        qSlider.setValue(7.0);
        add(qSlider);



        JPanel fsPanel = (JPanel) plotFs.graphic();
        fsPanel.setLayout(new GridBagLayout());
        GridBagConstraints gbcFs = new GridBagConstraints();
        gbcFs.gridx = 0;
        gbcFs.gridy = 0;
        gbcFs.gridheight = 1;
        fsPanel.add(plotFs.getPlot(), gbcFs);
        gbcFs.insets = new Insets(20, 0, 0, 0);
        gbcFs.gridheight = 1;
        gbcFs.gridx = 1;
        fsPanel.add(qSlider.graphic(), gbcFs);
        gbcFs.gridy = 1;
        gbcFs.gridx = gbcFs.gridy = 0;
        gbcFs.gridheight = 1;
        gbcFs.insets = new Insets(0, 0, 0, 0);






        meterMSD.addMSDSink(dsMSDcorP);
        DisplayPlot plotMSDcorUP = new DisplayPlot();
        plotMSDcorUP.setLabel("MSD cor P");
        plotMSDcorUP.getPlot().setXLog(true);
        plotMSDcorUP.setDoLegend(false);
        DataPumpListener pumpMSDcorP = new DataPumpListener(dsMSDcorP, plotMSDcorUP.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpMSDcorP);

        if (dsMSDcorU != null) {
            meterMSD.addMSDSink(dsMSDcorU);
            plotMSDcorUP.setLegend(new DataTag[]{dsMSDcorU.getTag()}, "U");
            DataPumpListener pumpMSDcorU = new DataPumpListener(dsMSDcorU, plotMSDcorUP.getDataSet().makeDataSink(), 1000);
            sim.integrator.getEventManager().addListener(pumpMSDcorU);
        }

        meterMSD.addMSDSink(dsPMSDhistory);
        DisplayPlot plotPMSDscatter = new DisplayPlot();
        plotPMSDscatter.setLabel("MSD vs. P");
        plotPMSDscatter.setXLabel("P");
        DataPumpListener pumpPMSDhistory = new DataPumpListener(dsPMSDhistory, plotPMSDscatter.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpPMSDhistory);
        plotPMSDscatter.setDoDrawLines(new DataTag[]{dsPMSDhistory.getTag()}, false);
        plotPMSDscatter.getPlot().setYLabel("MSD");
        plotPMSDscatter.setDoLegend(false);


        DisplayPlot plotHistogramP = new DisplayPlot();
        DataPumpListener pumpHistogramP = new DataPumpListener(dsHistogramP, plotHistogramP.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpHistogramP);
        plotHistogramP.setXLabel("P");
        plotHistogramP.setDoLegend(false);
        plotHistogramP.getPlot().setYLog(true);

        DataSourceHistogramMSD dsHistogramMSD = new DataSourceHistogramMSD(sim.integrator);
        meterMSD.addMSDSink(dsHistogramMSD);
        DisplayPlot plotHistogramMSD = new DisplayPlot();
        DataPumpListener pumpHistogramMSD = new DataPumpListener(dsHistogramMSD, plotHistogramMSD.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpHistogramMSD);
        plotHistogramMSD.setXLabel("MSD");
        plotHistogramMSD.setDoLegend(false);
        plotHistogramMSD.getPlot().setYLog(true);

        AccumulatorAverageFixed accSFac = new AccumulatorAverageFixed(1);  // just average, no uncertainty
        DataClusterer sfacClusterer = new DataClusterer(100, sim.getRandom());
        DeviceButtonGroup sfacButtons = null;
        DisplayPlot plotSFac = null;
        MeterStructureFactor meterSFacCluster = null;
        if (sim.box.getLeafList().size() <= 500) {
            MeterStructureFactor meterSFac = new MeterStructureFactor(space, sim.box, 15);
            accSFac.setPushInterval(1);
            DataPumpListener pumpSFac = new DataPumpListener(meterSFac, accSFac, 1000);
            sim.integrator.getEventManager().addListener(pumpSFac);
            dataStreamPumps.add(pumpSFac);
            plotSFac = new DisplayPlot();
            accSFac.addDataSink(plotSFac.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{accSFac.AVERAGE});
            plotSFac.setLabel("SFac");
            plotSFac.setDoDrawLines(new DataTag[]{meterSFac.getTag()}, false);
            meterSFac.setAtomTypeFactor(sim.speciesB.getAtomType(0), -1);

            double L = sim.box.getBoundary().getBoxSize().getX(0);
            double cut = 2.0 * Math.PI / L;

            meterSFacCluster = new MeterStructureFactor(space, sim.box, 3 * cut + 0.001);
            DataPumpListener pumpSFacCluster = new DataPumpListener(meterSFacCluster, sfacClusterer, 10);
            sim.integrator.getEventManager().addListener(pumpSFacCluster);
            pFork.addDataSink(sfacClusterer.makePressureSink());

            sfacButtons = new DeviceButtonGroup(sim.getController(), 5);
            sfacButtons.setLabel("B signal");
            AtomType typeB = sim.speciesB.getLeafType();
            MeterStructureFactor[] meters = new MeterStructureFactor[]{meterSFac, meterSFacCluster};
            sfacButtons.addButton("+1", new SFacButtonAction(meters, accSFac, sfacClusterer, typeB, +1));
            sfacButtons.addButton("-1", new SFacButtonAction(meters, accSFac, sfacClusterer, typeB, -1));
            double vB = sim.getSpace().powerD(sim.sigmaB);
            sfacButtons.addButton("+v", new SFacButtonAction(meters, accSFac, sfacClusterer, typeB, vB));
            sfacButtons.addButton("-v", new SFacButtonAction(meters, accSFac, sfacClusterer, typeB, -vB));
            sfacButtons.addButton("0", new SFacButtonAction(meters, accSFac, sfacClusterer, typeB, 0));
            sfacButtons.setSelected("+v");
        }

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
                    configStorageMSDPerc.reset();
                    configStorageMSDPerc.setEnabled(false);
                    meterPerc.reset();
                    dsMSDcorP.setEnabled(false);
                    dsHistogramP.setEnabled(false);
                    dsPMSDhistory.setEnabled(false);
                    dpAutocor.reset();
                    dpxyAutocor.reset();
                    pTensorAccum.setEnabled(false);
                    pTensorAccum.reset();
                    if (dsMSDcorU != null) dsMSDcorU.setEnabled(false);
                    meterFs.reset();
                    meterFsA.reset();
                    meterFsB.reset();
                    meterFOld.reset();
                    meterF.reset();
                    meterCorrelationAA.reset();
                    meterCorrelationAB.reset();
                    meterCorrelationBB.reset();
                    meterCorrelation.reset();
                    meterCorrelationPerp.reset();
                    meterCorrelationPar.reset();
                    meterCorrelationMag.reset();
                    diameterHash.setFac(1.0);
                    accSFac.reset();
                    sfacClusterer.reset();
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
                    configStorageMSDPerc.reset();
                    configStorageMSDPerc.setEnabled(true);
                    meterPerc.reset();
                    dsMSDcorP.setEnabled(true);
                    dsHistogramP.setEnabled(true);
                    dsPMSDhistory.setEnabled(true);
                    if (dsMSDcorU != null) dsMSDcorU.setEnabled(true);
                    if (colorCheckbox.getState()) dbox.setColorScheme(colorSchemeDeviation);
                    else if (colorDirectionCheckbox.getState()) dbox.setColorScheme(colorSchemeDirection);
                    if (showDispCheckbox.getState()) canvas.setDrawDisplacement(true);
                    if (flipDispCheckbox.getState()) canvas.setFlipDisplacement(true);
                    dpAutocor.reset();
                    dpxyAutocor.reset();
                    pTensorAccum.reset();
                    pTensorAccum.setEnabled(true);
                    meterFs.reset();
                    meterFsA.reset();
                    meterFsB.reset();
                    meterF.reset();
                    meterFOld.reset();
                    meterCorrelationAA.reset();
                    meterCorrelationAB.reset();
                    meterCorrelationBB.reset();
                    meterCorrelation.reset();
                    meterCorrelationPerp.reset();
                    meterCorrelationPar.reset();
                    meterCorrelationMag.reset();
                    diameterHash.setFac(showDispCheckbox.getState() ? 0.5 : 1.0);
                    accSFac.reset();
                    sfacClusterer.reset();
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

        DeviceSlider temperatureSelect = null;
        if (sim.potentialChoice != SimGlass.PotentialChoice.HS) {
            temperatureSelect = new DeviceSlider(sim.getController(), sim.integrator, "temperature");
            temperatureSelect.setPrecision(3);
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
        }
        configStorage.reset();
        configStorageLinear.reset();
        configStorageMSD.reset();
        configStorageMSDPerc.reset();
        dbox.repaint();

        IAction resetAction = new IAction() {
            public void actionPerformed() {
                rdfMeter.reset();

                getDisplayBox(sim.box).graphic().repaint();
            }
        };

        this.getController().getReinitButton().setPostAction(resetAction);
        this.getController().getResetAveragesButton().setPostAction(resetAction);

        if (temperatureSelect != null) getPanel().controlPanel.add(temperatureSelect.graphic(), vertGBC);

        add(rdfPlot);
        if (sim.potentialChoice != SimGlass.PotentialChoice.HS) {add(ePlot); add(tPlot);}
        add(densityBox);
        add(pDisplay);
        if (sim.potentialChoice != SimGlass.PotentialChoice.HS) {
            add(peDisplay);
        }
        JPanel plotPMSDPanel = new JPanel();
        plotPMSDPanel.setLayout(new BoxLayout(plotPMSDPanel, BoxLayout.Y_AXIS));
        plotPMSDPanel.add(plotMSDcorUP.graphic());
        DeviceSlider sliderPMSDinterval = new DeviceSlider(sim.getController(), new Modifier() {
            @Override
            public void setValue(double newValue) {
                dsHistogramMSD.setInterval((int) newValue);
                dsHistogramP.setInterval((int) newValue);
                dsPMSDhistory.setInterval((int) newValue);
                pumpHistogramMSD.actionPerformed();
                pumpHistogramP.actionPerformed();
                pumpPMSDhistory.actionPerformed();
            }

            @Override
            public double getValue() {
                return (double) dsHistogramP.getInterval();
            }

            @Override
            public Dimension getDimension() {
                return Null.DIMENSION;
            }

            @Override
            public String getLabel() {
                return "log2(interval)";
            }
        });
        sliderPMSDinterval.setMinimum(0);
        sliderPMSDinterval.setMaximum(30);
        sliderPMSDinterval.setNMajor(5);
        sliderPMSDinterval.setShowBorder(true);
        sliderPMSDinterval.setLabel("log2(histogram interval)");
        plotPMSDPanel.add(sliderPMSDinterval.graphic());
        plotPMSDPanel.add(plotHistogramP.graphic());
        plotPMSDPanel.add(plotHistogramMSD.graphic());
        plotPMSDPanel.add(plotPMSDscatter.graphic());
        JScrollPane plotsPane = new JScrollPane(plotPMSDPanel);
        java.awt.Dimension d = plotHistogramMSD.getPlot().getPreferredSize();
        d.height = 600;
        plotsPane.setPreferredSize(d);
        getPanel().tabbedPane.add("P cor MSD", plotsPane);

        if (sfacButtons != null) {
            JPanel plotSFacPanel = new JPanel();
            plotSFacPanel.setLayout(new BoxLayout(plotSFacPanel, BoxLayout.Y_AXIS));
            plotSFacPanel.add(sfacButtons.graphic());
            plotSFacPanel.add(plotSFac.graphic());
            JPanel sfacWidgetPanel = new JPanel();

            DeviceBox clusterIterBox = new DeviceBox();
            clusterIterBox.setController(sim.getController());
            clusterIterBox.setModifier(new ModifierGeneral(sfacClusterer, "maxIterations"));
            clusterIterBox.setInteger(true);
            clusterIterBox.setPrecision(0);
            clusterIterBox.setLabel("iterations");
            sfacWidgetPanel.add(clusterIterBox.graphic());
            DeviceButton clusterButton = new DeviceButton(sim.getController(), new IAction() {
                @Override
                public void actionPerformed() {
                    sfacClusterer.findClusters();
                }
            });
            clusterButton.setLabel("Cluster");
            sfacWidgetPanel.add(clusterButton.graphic());
            DeviceSlider cutClusterSlider = new DeviceSlider(sim.getController(), meterSFacCluster, "cutoff");
            cutClusterSlider.setLabel("SFac cutoff");
            cutClusterSlider.setShowBorder(true);
            cutClusterSlider.setMinimum(0);
            cutClusterSlider.setPrecision(2);
            cutClusterSlider.setMaximum(15);
            cutClusterSlider.setNMajor(5);
            cutClusterSlider.setShowValues(true);
            cutClusterSlider.setEditValues(true);
            sfacWidgetPanel.add(cutClusterSlider.graphic());
            plotSFacPanel.add(sfacWidgetPanel);
            sfacWidgetPanel = new JPanel();
            DeviceSlider nbrClusterSlider = new DeviceSlider(sim.getController(), sfacClusterer, "clusterNeighborDistance");
            nbrClusterSlider.setLabel("Cluster nbr distance");
            nbrClusterSlider.setShowBorder(true);
            nbrClusterSlider.setMinimum(0);
            nbrClusterSlider.setPrecision(2);
            nbrClusterSlider.setMaximum(3);
            nbrClusterSlider.setNMajor(3);
            nbrClusterSlider.setShowValues(true);
            nbrClusterSlider.setEditValues(true);
            sfacWidgetPanel.add(nbrClusterSlider.graphic());
            DeviceSlider nClusterSlider = new DeviceSlider(sim.getController(), sfacClusterer, "numClusters");
            nClusterSlider.setLabel("# of clusters");
            nClusterSlider.setShowBorder(true);
            nClusterSlider.setMinimum(0);
            nClusterSlider.setMaximum(5000);
            nClusterSlider.setPrecision(0);
            nClusterSlider.setNMajor(5);
            nClusterSlider.setShowValues(true);
            nClusterSlider.setEditValues(true);
            sfacWidgetPanel.add(nClusterSlider.graphic());
            plotSFacPanel.add(sfacWidgetPanel);

            JScrollPane plotsPaneSFac = new JScrollPane(plotSFacPanel);
            d = plotSFac.getPlot().getPreferredSize();
            d.height = 600;
            plotsPaneSFac.setPreferredSize(d);
            getPanel().tabbedPane.add("SFac", plotsPaneSFac);
        }

    }

    public static void main(String[] args) {
        SimGlass.GlassParams params = new SimGlass.GlassParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.doSwap = true;
            params.potential = SimGlass.PotentialChoice.HS;
            params.nA = 100;
            params.nB = 100;
            params.density = 1.6;
            params.D = 3;
        }
        SimGlass sim = new SimGlass(params.D, params.nA, params.nB, params.density, params.temperature, params.doSwap, params.potential, params.tStep);

        GlassGraphic ljmdGraphic = new GlassGraphic(sim);
        SimulationGraphic.makeAndDisplayFrame
                (ljmdGraphic.getPanel(), sim.potentialChoice+APP_NAME);
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

    public static class SFacButtonAction implements IAction {
        protected final double value;
        protected final AtomType type;
        protected final MeterStructureFactor[] meters;
        protected final AccumulatorAverageFixed acc;
        protected final DataClusterer sfacClusterer;

        public SFacButtonAction(MeterStructureFactor[] meters, AccumulatorAverageFixed acc, DataClusterer sfacClusterer, AtomType type, double value) {
            this.meters = meters;
            this.acc = acc;
            this.sfacClusterer = sfacClusterer;
            this.type = type;
            this.value = value;
        }

        public void actionPerformed() {
            for (MeterStructureFactor m : meters) {
                m.setAtomTypeFactor(type, value);
            }
            acc.reset();
            sfacClusterer.reset();
        }
    }
}