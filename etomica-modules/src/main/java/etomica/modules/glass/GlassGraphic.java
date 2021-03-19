/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.glass;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
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
import etomica.space.Vector;
import etomica.units.SimpleUnit;
import etomica.units.Unit;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;
import etomica.util.ParseArgs;
import net.miginfocom.swing.MigLayout;
import org.knowm.xchart.XYSeries;

import javax.swing.*;
import java.awt.*;
import java.awt.event.WindowEvent;
import java.util.ArrayList;
import java.util.List;

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
        DisplayPlotXChart rdfPlot = new DisplayPlotXChart();
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

        ConfigurationStorage configStorage = new ConfigurationStorage(sim.box, ConfigurationStorage.StorageType.LOG2);
        configStorage.setEnabled(false); // start isothermal
        sim.integrator.getEventManager().addListener(configStorage);
        DisplayBox dbox;
        DisplayBoxCanvasGlass canvas;
        DisplayCanvas c;
        dbox = new DisplayBox(sim.getController(), sim.box);
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
        AtomTestDeviation atomFilterDeviationPerc = new AtomTestDeviation(sim.box, configStorage);

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
                } else {
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

        DeviceSlider corIntervalSlider = new DeviceSlider(sim.getController(), null);
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
        MeterDensity densityMeter = new MeterDensity(sim.box);
        final DisplayTextBox densityBox = new DisplayTextBox();
        final DataPump densityPump = new DataPump(densityMeter, densityBox);
        IntegratorListenerAction densityPumpListener = new IntegratorListenerAction(densityPump);
        sim.integrator.getEventManager().addListener(densityPumpListener);
        densityPumpListener.setInterval(10);
        dataStreamPumps.add(densityPump);
        densityBox.setLabel("Number Density");

        AccumulatorPTensor pTensorAccum = new AccumulatorPTensor(sim.box, sim.integrator.getTimeStep());

        AccumulatorHistory energyHistory, peHistory;
        final AccumulatorAverageCollapsing peAccumulator;
        DisplayPlotXChart ePlot, tPlot;
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
            sim.integrator.getEventManager().addListener(kePump);
            keHistory.setPushInterval(5);
            dataStreamPumps.add(kePump);

            ePlot = new DisplayPlotXChart();
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

            tPlot = new DisplayPlotXChart();
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
        if (sim.box.getLeafList().size() > 200 && sim.potentialChoice != SimGlass.PotentialChoice.HS && false) {
            pTensorFork.addDataSink(dpAutocor);
        }

        //normalized AC of shear stress components
        AccumulatorAutocorrelationShearStress dpxyAutocor = new AccumulatorAutocorrelationShearStress(256, sim.integrator.getTimeStep());
        if (sim.box.getLeafList().size() > 200 && sim.potentialChoice != SimGlass.PotentialChoice.HS && false) {
            pTensorFork.addDataSink(dpxyAutocor);
        }


        DisplayPlotXChart plotPTensorAccum = new DisplayPlotXChart();
        plotPTensorAccum.setLabel("viscosity(t)");
        plotPTensorAccum.setDoLegend(false);
        plotPTensorAccum.setXLog(true);
        plotPTensorAccum.setYLog(true);
        add(plotPTensorAccum);
        DataPumpListener pTensorAccumPump = new DataPumpListener(pTensorAccum, plotPTensorAccum.makeSink("viscosity"), 1000);
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
        DisplayPlotXChart plotP = new DisplayPlotXChart();
        historyP.addDataSink(plotP.getDataSet().makeDataSink());
        plotP.setLabel("P");
        plotP.setDoLegend(false);
        add(plotP);
        DataProcessorErrorBar pAutoCorErr = new DataProcessorErrorBar("err+");
        dpAutocor.getAvgErrFork().addDataSink(pAutoCorErr);

        DataSourceMSDcorP dsMSDcorP = new DataSourceMSDcorP(sim.timeSource);
        pFork.addDataSink(dsMSDcorP);
        DataSourceBlockAvgCor dsCorP = new DataSourceBlockAvgCor(sim.timeSource);
        pFork.addDataSink(dsCorP);

        DataSourceBlockAvgCor dsCorKE = new DataSourceBlockAvgCor(sim.timeSource);
        if (sim.potentialChoice != SimGlass.PotentialChoice.HS) {
            MeterKineticEnergy keMeter = new MeterKineticEnergy(sim.box);
            DataPumpListener kePump = new DataPumpListener(keMeter, dsCorKE, 1 << 3);
            sim.integrator.getEventManager().addListener(kePump);
            dsCorKE.setMinInterval(3);
        }
        DataSourceCorMSD dsCorMSD = new DataSourceCorMSD(sim.timeSource);
        dsCorMSD.setMinInterval(3);

        DataSourceHisogram dsHistogramP = new DataSourceHisogram(sim.timeSource);
        pFork.addDataSink(dsHistogramP);

        DataSourcePMSDHistory dsPMSDhistory = new DataSourcePMSDHistory(sim.timeSource);
        pFork.addDataSink(dsPMSDhistory);

        DataSourceMSDcorP dsMSDcorU;
        if (sim.potentialChoice != SimGlass.PotentialChoice.HS && false) {
            // disabled.  U is not differently correlated than P
            dsMSDcorU = new DataSourceMSDcorP(sim.timeSource, log2peInterval);
            peFork.addDataSink(dsMSDcorU);
        } else {
            dsMSDcorU = null;
        }



        //Gs: total
        int gsUpdateInterval = sim.getSpace().D() == 2 ? 10000 : 2000;
        double xGsMax = 3;
        int gsMinConfig = 6;
        MeterGs meterGs = new MeterGs(configStorage);
        meterGs.setMinConfigIndex(gsMinConfig);
        meterGs.setConfigIndex(12);
        configStorage.addListener(meterGs);
        meterGs.getXDataSource().setXMax(xGsMax);
        meterGs.getXDataSource().setNValues(150);
        DisplayPlotXChart gsPlot = new DisplayPlotXChart();
        DataPumpListener pumpGs = new DataPumpListener(meterGs, gsPlot.getDataSet().makeDataSink(), gsUpdateInterval);
        sim.integrator.getEventManager().addListener(pumpGs);
        gsPlot.setLabel(" Gs ");
        //Gs: A
        MeterGs meterGsA = new MeterGs(configStorage);
        meterGsA.setMinConfigIndex(gsMinConfig);
        meterGsA.setConfigIndex(12);
        meterGsA.setAtomTypes(sim.speciesA.getLeafType());
        configStorage.addListener(meterGsA);
        meterGsA.getXDataSource().setXMax(xGsMax);
        meterGsA.getXDataSource().setNValues(150);
        DataPumpListener pumpGsA = new DataPumpListener(meterGsA, gsPlot.getDataSet().makeDataSink(), gsUpdateInterval);
        sim.integrator.getEventManager().addListener(pumpGsA);

        //Gs: B
        MeterGs meterGsB = new MeterGs(configStorage);
        meterGsB.setMinConfigIndex(gsMinConfig);
        meterGsB.setConfigIndex(12);
        meterGsB.setAtomTypes(sim.speciesB.getLeafType());
        configStorage.addListener(meterGsB);
        meterGsB.getXDataSource().setXMax(xGsMax);
        meterGsB.getXDataSource().setNValues(150);
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
                meterGs.setConfigIndex(log2prevConfig);
                meterGsA.setConfigIndex(log2prevConfig);
                meterGsB.setConfigIndex(log2prevConfig);
                pumpGs.actionPerformed();
                pumpGsA.actionPerformed();
                pumpGsB.actionPerformed();
            }

            @Override
            public double getValue() {
                return meterGs.getConfigIndex();
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
        gsPrevSampleSlider.setMaximum(30);
        gsPrevSampleSlider.setMinimum(0);
        gsPrevSampleSlider.setLabel("log2(previous sample)");

        JPanel gsPanel = new JPanel();
        gsPanel.setLayout(new MigLayout("fill"));
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridheight = 1;
        gsPanel.add(gsPlot.getPanel(), "wrap, grow");
        gbc.insets = new Insets(20, 0, 0, 0);
        gbc.gridheight = 1;
        gbc.gridx = 1;
        gsPanel.add(gsPrevSampleSlider.graphic(), "alignx center");
        gbc.gridy = 1;
        gbc.gridx = gbc.gridy = 0;
        gbc.gridheight = 1;
        gbc.insets = new Insets(0, 0, 0, 0);
        this.addAsTab(gsPanel, gsPlot.getLabel(), true);


        //Strings
        DataSourceStrings meterStrings = new DataSourceStrings(configStorage, 3);
        configStorage.addListener(meterStrings);
        DisplayPlotXChart plotStrings = new DisplayPlotXChart();
        DataPumpListener pumpStrings = new DataPumpListener(meterStrings, plotStrings.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpStrings);
        plotStrings.setLabel("strings");
        plotStrings.getPlot().setXLog(true);
        plotStrings.setDoLegend(false);
        add(plotStrings);

        //Percolation
        atomFilterDeviationPerc.setDoMobileOnly(false);
        DataSourcePercolation meterPerc = new DataSourcePercolation(configStorage, atomFilterDeviationPerc, 8);
        configStorage.addListener(meterPerc);
        DisplayPlotXChart plotPerc = new DisplayPlotXChart();
        DataPumpListener pumpPerc = new DataPumpListener(meterPerc, plotPerc.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpPerc);

        DataSourcePercolation.ImmFractionSource meterImmFrac = meterPerc.makeImmFractionSource();
        DataPumpListener pumpImmFrac = new DataPumpListener(meterImmFrac, plotPerc.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpImmFrac);
        DataSourcePercolation.ImmFractionByTypeSource meterImmFracA = meterPerc.makeImmFractionSource(sim.speciesA.getLeafType());
        DataPumpListener pumpImmFracA = new DataPumpListener(meterImmFracA, plotPerc.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpImmFracA);
        DataSourcePercolation.ImmFractionByTypeSource meterImmFracB = meterPerc.makeImmFractionSource(sim.speciesB.getLeafType());
        DataPumpListener pumpImmFracB = new DataPumpListener(meterImmFracB, plotPerc.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpImmFracB);

        DataSourceQ4 meterQ4 = new DataSourceQ4(configStorage, 8);
        configStorage.addListener(meterQ4);
        DataPumpListener pumpQ4 = new DataPumpListener(meterQ4, plotPerc.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpQ4);

        plotPerc.setLegend(new DataTag[]{meterPerc.getTag()}, "perc. prob.");
        plotPerc.setLegend(new DataTag[]{meterImmFrac.getTag()}, "imm. frac.");
        plotPerc.setLegend(new DataTag[]{meterImmFracA.getTag()}, "imm. frac. A");
        plotPerc.setLegend(new DataTag[]{meterImmFracB.getTag()}, "imm. frac. B");
        plotPerc.setLegend(new DataTag[]{meterQ4.getTag()}, "Q4");
        plotPerc.setLabel("perc");
        plotPerc.getPlot().setXLog(true);

        DisplayPlotXChart plotChi4 = new DisplayPlotXChart();
        plotChi4.getPlot().setXLog(true);
        DataSourcePercolation.Chi4Source meterChi4Star = meterPerc.makeChi4Source();
        DataPumpListener pumpChi4Star = new DataPumpListener(meterChi4Star, plotChi4.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpChi4Star);
        plotChi4.setLegend(new DataTag[]{meterChi4Star.getTag()}, "chi4*");
        DataSourceQ4.MeterChi4 meterChi4 = meterQ4.makeChi4Meter();
        DataPumpListener pumpChi4 = new DataPumpListener(meterChi4, plotChi4.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpChi4);
        plotChi4.setLegend(new DataTag[]{meterChi4.getTag()}, "chi4");

        JPanel plotPercPanel = new JPanel(new MigLayout("fill"));
        JScrollPane plotsPerc = new JScrollPane(plotPercPanel);
        getPanel().tabbedPane.add("perc", plotsPerc);


        //Percolation slider
        DeviceSlider percDrSlider = new DeviceSlider(sim.getController(), new Modifier() {
            @Override
            public void setValue(double newValue) {
                if (newValue == getValue()) return;
                atomFilterDeviationPerc.setMinDistance(newValue);
                meterPerc.zeroData();
                meterQ4.setMaxDr(newValue);
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

        plotPercPanel.add(plotPerc.getPanel());
        plotPercPanel.add(plotChi4.getPanel(), "wrap");
        plotPercPanel.add(percDrSlider.graphic(), "alignx center");
        plotPercPanel.add(percMinLog2StepSlider.getPanel(), "alignx center");

        DataSourcePercolation0 meterPerc0 = new DataSourcePercolation0(sim.box, sim.getRandom());
        meterPerc0.setImmFracs(new double[]{0.05, 0.1, 0.15, 0.20, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.8, 0.85, 0.9, 0.95, 1});
        AccumulatorAverageFixed accPerc0 = new AccumulatorAverageFixed(10);
        accPerc0.setPushInterval(1);
        DataPumpListener pumpPerc0 = new DataPumpListener(meterPerc0, accPerc0, 1000);
        DisplayPlotXChart plotPerc0 = new DisplayPlotXChart();
        plotPerc0.setXLabel("Atom Fraction");
        plotPerc0.setYLabel("Percolation Fraction");
        accPerc0.addDataSink(plotPerc0.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{accPerc0.AVERAGE});
        plotPerc0.setLegend(new DataTag[]{meterPerc0.getTag()}, "random");
        sim.integrator.getEventManager().addListener(pumpPerc0);
        dataStreamPumps.add(pumpPerc0);
        plotPerc0.setLabel("perc(imm)");
        add(plotPerc0);
        DataSourcePercolation.PercolationByImmFrac meterImmFracPerc = meterPerc.makePerclationByImmFracSource();
        DataPumpListener pumpImmFracPerc = new DataPumpListener(meterImmFracPerc, plotPerc0.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpImmFracPerc);
        plotPerc0.setLegend(new DataTag[]{meterImmFracPerc.getTag()}, "immobile");

        DisplayPlotXChart plotCorSelf = new DisplayPlotXChart();
        plotCorSelf.getPlot().setXLog(true);
        plotCorSelf.setDoLegend(true);

        MeterCorrelationSelf meterCorrelationSelf = new MeterCorrelationSelf(configStorage);
        configStorage.addListener(meterCorrelationSelf);
        DataPumpListener pumpCorSelf = new DataPumpListener(meterCorrelationSelf, plotCorSelf.getDataSet().makeDataSink(), 1000);
        plotCorSelf.setLegend(new DataTag[]{meterCorrelationSelf.getTag()}, "dot");
        sim.integrator.getEventManager().addListener(pumpCorSelf);
        MeterCorrelationSelf meterCorrelationSelfA = new MeterCorrelationSelf(configStorage);
        meterCorrelationSelfA.setAtomType(sim.speciesA.getLeafType());
        configStorage.addListener(meterCorrelationSelfA);
        DataPumpListener pumpCorSelfA = new DataPumpListener(meterCorrelationSelfA, plotCorSelf.getDataSet().makeDataSink(), 1000);
        plotCorSelf.setLegend(new DataTag[]{meterCorrelationSelfA.getTag()}, "dotA");
        sim.integrator.getEventManager().addListener(pumpCorSelfA);
        MeterCorrelationSelf meterCorrelationSelfB = new MeterCorrelationSelf(configStorage);
        meterCorrelationSelfB.setAtomType(sim.speciesB.getLeafType());
        configStorage.addListener(meterCorrelationSelfB);
        DataPumpListener pumpCorSelfB = new DataPumpListener(meterCorrelationSelfB, plotCorSelf.getDataSet().makeDataSink(), 1000);
        plotCorSelf.setLegend(new DataTag[]{meterCorrelationSelfB.getTag()}, "dotB");
        sim.integrator.getEventManager().addListener(pumpCorSelfB);

        MeterCorrelationSelf meterCorrelationSelfMagA = new MeterCorrelationSelf(configStorage, MeterCorrelationSelf.CorrelationType.MAGNITUDE);
        meterCorrelationSelfMagA.setAtomType(sim.speciesA.getLeafType());
        configStorage.addListener(meterCorrelationSelfMagA);
        DataPumpListener pumpCorSelfMagA = new DataPumpListener(meterCorrelationSelfMagA, plotCorSelf.getDataSet().makeDataSink(), 1000);
        plotCorSelf.setLegend(new DataTag[]{meterCorrelationSelfMagA.getTag()}, "|r|A");
        sim.integrator.getEventManager().addListener(pumpCorSelfMagA);
        MeterCorrelationSelf meterCorrelationSelfMagB = new MeterCorrelationSelf(configStorage, MeterCorrelationSelf.CorrelationType.MAGNITUDE);
        meterCorrelationSelfMagB.setAtomType(sim.speciesB.getLeafType());
        configStorage.addListener(meterCorrelationSelfMagB);
        DataPumpListener pumpCorSelfMagB = new DataPumpListener(meterCorrelationSelfMagB, plotCorSelf.getDataSet().makeDataSink(), 1000);
        plotCorSelf.setLegend(new DataTag[]{meterCorrelationSelfMagB.getTag()}, "|r|B");
        sim.integrator.getEventManager().addListener(pumpCorSelfMagB);

        CorrelationSelf2 correlationSelf2 = new CorrelationSelf2(configStorage, CorrelationSelf2.CorrelationType.TOTAL, 0.001, 20);
        configStorage.addListener(correlationSelf2);
        DisplayPlotXChart plotCorSelf2 = new DisplayPlotXChart();
        for (int i = 2; i < correlationSelf2.getNumDt(); i++) {
            CorrelationSelf2.MeterCorrelationSelf2 m = correlationSelf2.makeMeter(i);
            DataPumpListener p = new DataPumpListener(m, plotCorSelf2.getDataSet().makeDataSink(), 1000);
            plotCorSelf2.setLegend(new DataTag[]{m.getTag()}, "" + i);
            sim.integrator.getEventManager().addListener(p);
        }

        JPanel plotCorSelfPanel = new JPanel(new MigLayout("fill"));
        plotCorSelfPanel.add(plotCorSelf.graphic(), "grow");
        plotCorSelfPanel.add(plotCorSelf2.graphic(), "grow");
        JScrollPane plotsCorSelfPane = new JScrollPane(plotCorSelfPanel);
        getPanel().tabbedPane.add("cor self", plotsCorSelfPane);


        int corUpdateInterval = sim.getSpace().D() == 2 ? 10000 : 2000;
        double xCorMax = 5;

        int minCorSample = 7;
        MeterCorrelation meterCorrelation = new MeterCorrelation(configStorage);
        meterCorrelation.setMinPrevSample(minCorSample);
        meterCorrelation.setPrevSampleIndex(7);
        configStorage.addListener(meterCorrelation);
        meterCorrelation.getXDataSource().setXMax(xCorMax);
        DisplayPlotXChart correlationPlot = new DisplayPlotXChart();
        DataPumpListener pumpCorrelation = new DataPumpListener(meterCorrelation, correlationPlot.getDataSet().makeDataSink(), corUpdateInterval);
        sim.integrator.getEventManager().addListener(pumpCorrelation);
        correlationPlot.setLabel("cor");

        MeterCorrelation meterCorrelationAA = new MeterCorrelation(configStorage);
        meterCorrelationAA.setAtomTypes(sim.speciesA.getLeafType(), sim.speciesA.getLeafType());
        meterCorrelationAA.setPrevSampleIndex(7);
        meterCorrelationAA.setMinPrevSample(minCorSample);
        configStorage.addListener(meterCorrelationAA);
        meterCorrelationAA.getXDataSource().setXMax(xCorMax);
        DataPumpListener pumpCorrelationAA = new DataPumpListener(meterCorrelationAA, correlationPlot.getDataSet().makeDataSink(), corUpdateInterval);
        sim.integrator.getEventManager().addListener(pumpCorrelationAA);
        MeterCorrelation meterCorrelationAB = new MeterCorrelation(configStorage);
        meterCorrelationAB.setAtomTypes(sim.speciesA.getLeafType(), sim.speciesB.getLeafType());
        meterCorrelationAB.setPrevSampleIndex(7);
        meterCorrelationAB.setMinPrevSample(minCorSample);
        configStorage.addListener(meterCorrelationAB);
        meterCorrelationAB.getXDataSource().setXMax(xCorMax);
        DataPumpListener pumpCorrelationAB = new DataPumpListener(meterCorrelationAB, correlationPlot.getDataSet().makeDataSink(), corUpdateInterval);
        sim.integrator.getEventManager().addListener(pumpCorrelationAB);
        MeterCorrelation meterCorrelationBB = new MeterCorrelation(configStorage);
        meterCorrelationBB.setAtomTypes(sim.speciesB.getLeafType(), sim.speciesB.getLeafType());
        meterCorrelationBB.setPrevSampleIndex(7);
        meterCorrelationBB.setMinPrevSample(minCorSample);
        configStorage.addListener(meterCorrelationBB);
        meterCorrelationBB.getXDataSource().setXMax(xCorMax);
        DataPumpListener pumpCorrelationBB = new DataPumpListener(meterCorrelationBB, correlationPlot.getDataSet().makeDataSink(), corUpdateInterval);
        sim.integrator.getEventManager().addListener(pumpCorrelationBB);

        correlationPlot.setLegend(new DataTag[]{meterCorrelation.getTag()}, "total");
        correlationPlot.setLegend(new DataTag[]{meterCorrelationAA.getTag()}, "AA");
        correlationPlot.setLegend(new DataTag[]{meterCorrelationAB.getTag()}, "AB");
        correlationPlot.setLegend(new DataTag[]{meterCorrelationBB.getTag()}, "BB");

        MeterCorrelation meterCorrelationPerp = new MeterCorrelation(configStorage, MeterCorrelation.CorrelationType.PERPENDICULAR);
        meterCorrelationPerp.setPrevSampleIndex(7);
        meterCorrelationPerp.setMinPrevSample(minCorSample);
        configStorage.addListener(meterCorrelationPerp);
        meterCorrelationPerp.getXDataSource().setXMax(xCorMax);
        DataPumpListener pumpCorrelationPerp = new DataPumpListener(meterCorrelationPerp, correlationPlot.getDataSet().makeDataSink(), corUpdateInterval);
        sim.integrator.getEventManager().addListener(pumpCorrelationPerp);
        correlationPlot.setLegend(new DataTag[]{meterCorrelationPerp.getTag()}, "_|_");
        MeterCorrelation meterCorrelationPar = new MeterCorrelation(configStorage, MeterCorrelation.CorrelationType.PARALLEL);
        meterCorrelationPar.setPrevSampleIndex(7);
        meterCorrelationPar.setMinPrevSample(minCorSample);
        configStorage.addListener(meterCorrelationPar);
        meterCorrelationPar.getXDataSource().setXMax(xCorMax);
        DataPumpListener pumpCorrelationPar = new DataPumpListener(meterCorrelationPar, correlationPlot.getDataSet().makeDataSink(), corUpdateInterval);
        sim.integrator.getEventManager().addListener(pumpCorrelationPar);
        correlationPlot.setLegend(new DataTag[]{meterCorrelationPar.getTag()}, "||");

        MeterCorrelation meterCorrelationMag = new MeterCorrelation(configStorage, MeterCorrelation.CorrelationType.MAGNITUDE);
        meterCorrelationMag.setPrevSampleIndex(7);
        meterCorrelationMag.setMinPrevSample(minCorSample);
        configStorage.addListener(meterCorrelationMag);
        meterCorrelationMag.getXDataSource().setXMax(xCorMax);
        DataPumpListener pumpCorrelationMag = new DataPumpListener(meterCorrelationMag, correlationPlot.getDataSet().makeDataSink(), corUpdateInterval);
        sim.integrator.getEventManager().addListener(pumpCorrelationMag);
        correlationPlot.setLegend(new DataTag[]{meterCorrelationMag.getTag()}, "|r|");
        MeterCorrelation meterCorrelationAAMag = new MeterCorrelation(configStorage, MeterCorrelation.CorrelationType.MAGNITUDE);
        meterCorrelationAAMag.setAtomTypes(sim.speciesA.getLeafType(), sim.speciesA.getLeafType());
        meterCorrelationAAMag.setPrevSampleIndex(7);
        meterCorrelationAAMag.setMinPrevSample(minCorSample);
        configStorage.addListener(meterCorrelationAAMag);
        meterCorrelationAAMag.getXDataSource().setXMax(xCorMax);
        DataPumpListener pumpCorrelationAAMag = new DataPumpListener(meterCorrelationAAMag, correlationPlot.getDataSet().makeDataSink(), corUpdateInterval);
        sim.integrator.getEventManager().addListener(pumpCorrelationAAMag);
        correlationPlot.setLegend(new DataTag[]{meterCorrelationAAMag.getTag()}, "|r|AA");
        MeterCorrelation meterCorrelationABMag = new MeterCorrelation(configStorage, MeterCorrelation.CorrelationType.MAGNITUDE);
        meterCorrelationABMag.setAtomTypes(sim.speciesA.getLeafType(), sim.speciesB.getLeafType());
        meterCorrelationABMag.setPrevSampleIndex(7);
        meterCorrelationABMag.setMinPrevSample(minCorSample);
        configStorage.addListener(meterCorrelationABMag);
        meterCorrelationABMag.getXDataSource().setXMax(xCorMax);
        DataPumpListener pumpCorrelationABMag = new DataPumpListener(meterCorrelationABMag, correlationPlot.getDataSet().makeDataSink(), corUpdateInterval);
        sim.integrator.getEventManager().addListener(pumpCorrelationABMag);
        correlationPlot.setLegend(new DataTag[]{meterCorrelationABMag.getTag()}, "|r|AB");
        MeterCorrelation meterCorrelationBBMag = new MeterCorrelation(configStorage, MeterCorrelation.CorrelationType.MAGNITUDE);
        meterCorrelationBBMag.setAtomTypes(sim.speciesB.getLeafType(), sim.speciesB.getLeafType());
        meterCorrelationBBMag.setPrevSampleIndex(7);
        meterCorrelationBBMag.setMinPrevSample(minCorSample);
        configStorage.addListener(meterCorrelationBBMag);
        meterCorrelationBBMag.getXDataSource().setXMax(xCorMax);
        DataPumpListener pumpCorrelationBBMag = new DataPumpListener(meterCorrelationBBMag, correlationPlot.getDataSet().makeDataSink(), corUpdateInterval);
        sim.integrator.getEventManager().addListener(pumpCorrelationBBMag);
        correlationPlot.setLegend(new DataTag[]{meterCorrelationBBMag.getTag()}, "|r|BB");


        DeviceSlider corMinPrevSampleSlider = new DeviceSlider(sim.getController(), new Modifier() {
            @Override
            public void setValue(double newValue) {
                if (newValue == getValue()) return;
                int log2prevConfig = (int) Math.round(newValue);
                if (log2prevConfig > meterCorrelation.getPrevSampleIndex())
                    throw new IllegalArgumentException("min can't be greater than value");
                meterCorrelationAA.setMinPrevSample(log2prevConfig);
                meterCorrelationAB.setMinPrevSample(log2prevConfig);
                meterCorrelationBB.setMinPrevSample(log2prevConfig);
                meterCorrelation.setMinPrevSample(log2prevConfig);
                meterCorrelationPerp.setMinPrevSample(log2prevConfig);
                meterCorrelationPar.setMinPrevSample(log2prevConfig);
                meterCorrelationMag.setMinPrevSample(log2prevConfig);
                meterCorrelationAAMag.setMinPrevSample(log2prevConfig);
                meterCorrelationABMag.setMinPrevSample(log2prevConfig);
                meterCorrelationBBMag.setMinPrevSample(log2prevConfig);
                pumpCorrelation.actionPerformed();
                pumpCorrelationAA.actionPerformed();
                pumpCorrelationAB.actionPerformed();
                pumpCorrelationBB.actionPerformed();
                pumpCorrelationPar.actionPerformed();
                pumpCorrelationPerp.actionPerformed();
                pumpCorrelationMag.actionPerformed();
                pumpCorrelationAAMag.actionPerformed();
                pumpCorrelationABMag.actionPerformed();
                pumpCorrelationBBMag.actionPerformed();
            }

            @Override
            public double getValue() {
                return meterCorrelation.getMinPrevSample();
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
        corMinPrevSampleSlider.setShowBorder(true);
        corMinPrevSampleSlider.setShowValues(true);
        corMinPrevSampleSlider.setNMajor(5);
        corMinPrevSampleSlider.setMaximum(15);
        corMinPrevSampleSlider.setMinimum(0);
        corMinPrevSampleSlider.setLabel("log2(min previous sample)");

        DeviceSlider corPrevSampleSlider = new DeviceSlider(sim.getController(), new Modifier() {
            @Override
            public void setValue(double newValue) {
                if (newValue == getValue()) return;
                int log2prevConfig = (int) Math.round(newValue);
                if (log2prevConfig < meterCorrelation.getMinPrevSample())
                    throw new IllegalArgumentException("value can't be less than min");
                meterCorrelationAA.setPrevSampleIndex(log2prevConfig);
                meterCorrelationAB.setPrevSampleIndex(log2prevConfig);
                meterCorrelationBB.setPrevSampleIndex(log2prevConfig);
                meterCorrelation.setPrevSampleIndex(log2prevConfig);
                meterCorrelationPerp.setPrevSampleIndex(log2prevConfig);
                meterCorrelationPar.setPrevSampleIndex(log2prevConfig);
                meterCorrelationMag.setPrevSampleIndex(log2prevConfig);
                meterCorrelationAAMag.setPrevSampleIndex(log2prevConfig);
                meterCorrelationABMag.setPrevSampleIndex(log2prevConfig);
                meterCorrelationBBMag.setPrevSampleIndex(log2prevConfig);
                pumpCorrelation.actionPerformed();
                pumpCorrelationAA.actionPerformed();
                pumpCorrelationAB.actionPerformed();
                pumpCorrelationBB.actionPerformed();
                pumpCorrelationPar.actionPerformed();
                pumpCorrelationPerp.actionPerformed();
                pumpCorrelationMag.actionPerformed();
                pumpCorrelationAAMag.actionPerformed();
                pumpCorrelationABMag.actionPerformed();
                pumpCorrelationBBMag.actionPerformed();
            }

            @Override
            public double getValue() {
                return meterCorrelation.getPrevSampleIndex();
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
        corPrevSampleSlider.setMaximum(30);
        corPrevSampleSlider.setMinimum(0);
        corPrevSampleSlider.setLabel("log2(previous sample)");

        corIntervalSlider.setModifier(new Modifier() {
            @Override
            public void setValue(double newValue) {
                int log2interval = (int) Math.round(newValue);
                int interval = 1 << log2interval;
                configStorageLinear.setSampleInterval(interval);
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

        JPanel corPanel = new JPanel(new MigLayout("fill"));
        corPanel.add(correlationPlot.getPanel(), "span 1 2, grow");
        corPanel.add(corMinPrevSampleSlider.graphic(), "wrap");
        corPanel.add(corPrevSampleSlider.graphic(), "");
        addAsTab(corPanel, correlationPlot.getLabel(), true);


        //unnormalized AC of all stress tensor components
        DisplayPlotXChart plotPTensorAutocor = new DisplayPlotXChart();
        plotPTensorAutocor.setLabel("P Tensor autocor");
        plotPTensorAutocor.setLegend(new DataTag[]{dpAutocor.getTag()}, "avg");
        dpAutocor.addDataSink(plotPTensorAutocor.getDataSet().makeDataSink());
//        add(plotPTensorAutocor);
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
        JPanel ptacPanel = new JPanel();
        ptacPanel.setLayout(new GridBagLayout());
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridheight = 2;
        ptacPanel.add(plotPTensorAutocor.getPanel(), gbc);
        gbc.insets = new Insets(20, 0, 0, 0);
        gbc.gridheight = 1;
        gbc.gridx = 1;
        ptacPanel.add(nMaxSlider.getPanel(), gbc);
        gbc.gridy = 1;
        ptacPanel.add(pushIntervalSlider.getPanel(), gbc);

        pAutoCorErr.setDataSink(plotPTensorAutocor.getDataSet().makeDataSink());
        plotPTensorAutocor.setLegend(new DataTag[]{dpAutocor.getAvgErrFork().getTag(), pAutoCorErr.getTag()}, "err+");

        //normalized AC of shear stress components
        DisplayPlotXChart plotPxyTensorAutocor = new DisplayPlotXChart();
        plotPxyTensorAutocor.setLabel("Pxy autocor");
        plotPxyTensorAutocor.setLegend(new DataTag[]{dpxyAutocor.getTag()}, "avg");
        dpxyAutocor.addDataSink(plotPxyTensorAutocor.getDataSet().makeDataSink());
//        add(plotPxyTensorAutocor);
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
        JPanel shearacPanel = new JPanel();
        shearacPanel.setLayout(new GridBagLayout());
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridheight = 2;
        shearacPanel.add(plotPxyTensorAutocor.getPanel(), gbc);
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
            DisplayPlotXChart plotPE = new DisplayPlotXChart();
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
        DataSourceMSD meterMSD = new DataSourceMSD(configStorage);
        configStorage.addListener(meterMSD);
        DisplayPlotXChart plotMSD = new DisplayPlotXChart();
        DataPumpListener pumpMSD = new DataPumpListener(meterMSD, plotMSD.getDataSet().makeDataSink(), 1000);
        plotMSD.setLegend(new DataTag[]{meterMSD.getTag()}, "total");
        sim.integrator.getEventManager().addListener(pumpMSD);
        plotMSD.setLabel("MSD");
        plotMSD.getPlot().setYLog(true);
        plotMSD.getPlot().setXLog(true);
        add(plotMSD);
        DataSourceMSD meterMSDA = new DataSourceMSD(configStorage, sim.speciesA.getLeafType());
        configStorage.addListener(meterMSDA);
        DataPumpListener pumpMSDA = new DataPumpListener(meterMSDA, plotMSD.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpMSDA);
        plotMSD.setLegend(new DataTag[]{meterMSDA.getTag()}, "A");
        DataSourceMSD meterMSDB = new DataSourceMSD(configStorage, sim.speciesB.getLeafType());
        configStorage.addListener(meterMSDB);
        DataPumpListener pumpMSDB = new DataPumpListener(meterMSDB, plotMSD.getDataSet().makeDataSink(), 1000);
        plotMSD.setLegend(new DataTag[]{meterMSDB.getTag()}, "B");
        sim.integrator.getEventManager().addListener(pumpMSDB);

        //VAC
        configStorage.setDoVelocity(true);
        DataSourceVAC meterVAC = new DataSourceVAC(configStorage);
        configStorage.addListener(meterVAC);
        DisplayPlotXChart plotVAC = new DisplayPlotXChart();
        DataPumpListener pumpVAC = new DataPumpListener(meterVAC, plotVAC.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpVAC);
        plotVAC.setLabel("VAC");
        plotVAC.getPlot().setXLog(true);
        add(plotVAC);


        DataSourceAlpha2 meterAlpha2 = new DataSourceAlpha2(configStorage);
        configStorage.addListener(meterAlpha2);
        DisplayPlotXChart plotAlpha2 = new DisplayPlotXChart();
        DataPumpListener pumpAlpha2 = new DataPumpListener(meterAlpha2, plotAlpha2.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpAlpha2);
        plotAlpha2.setLabel("alpha2");
        plotAlpha2.getPlot().setXLog(true);
        plotAlpha2.setDoLegend(false);
        add(plotAlpha2);


        //F - new
        DataSourceF meterF = new DataSourceF(configStorage);
        configStorage.addListener(meterF);
        DisplayPlotXChart plotF = new DisplayPlotXChart();
        DataPumpListener pumpF = new DataPumpListener(meterF, plotF.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpF);
        plotF.setLabel("F");
        plotF.getPlot().setXLog(true);
        add(plotF);


        plotF.setLegend(new DataTag[]{meterF.getTag()}, "new");



        //Fs: TOTAL
        double qx = 7.0;
        Vector qVec = sim.getSpace().makeVector();
        qVec.setX(0, qx);

        DataSourceFs meterFs = new DataSourceFs(configStorage);
        meterFs.setQ(qVec);
        configStorage.addListener(meterFs);
        DisplayPlotXChart plotFs = new DisplayPlotXChart();
        DataPumpListener pumpFs = new DataPumpListener(meterFs, plotFs.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpFs);
        plotFs.setLabel("Fs");
        plotFs.getPlot().setXLog(true);

        //Fs: A
        DataSourceFs meterFsA = new DataSourceFs(configStorage);
        meterFsA.setQ(qVec);
        meterFsA.setAtomType(sim.speciesA.getLeafType());
        configStorage.addListener(meterFsA);
        DataPumpListener pumpFsA = new DataPumpListener(meterFsA, plotFs.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpFsA);

        //Fs: B
        DataSourceFs meterFsB = new DataSourceFs(configStorage);
        meterFsB.setQ(qVec);
        meterFsB.setAtomType(sim.speciesB.getLeafType());
        configStorage.addListener(meterFsB);
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


        JPanel fsPanel = new JPanel(new MigLayout("fill"));
        fsPanel.add(plotFs.getPanel(), "grow");
        fsPanel.add(qSlider.graphic(), "");
        this.addAsTab(fsPanel, plotFs.getLabel(), true);


        meterMSD.addMSDSink(dsMSDcorP);
        DisplayPlotXChart plotMSDcorUP = new DisplayPlotXChart();
        plotMSDcorUP.getPlot().setTitle("correlation");
        plotMSDcorUP.getPlot().setXLog(true);
        DataPumpListener pumpMSDcorP = new DataPumpListener(dsMSDcorP, plotMSDcorUP.getDataSet().makeDataSink(), 1000);
        plotMSDcorUP.setLegend(new DataTag[]{dsMSDcorP.getTag()}, "P-MSD");
        sim.integrator.getEventManager().addListener(pumpMSDcorP);
        meterMSD.addMSDSink(dsCorMSD);
        DataPumpListener pumpMSDcor = new DataPumpListener(dsCorMSD, plotMSDcorUP.getDataSet().makeDataSink(), 1000);
        plotMSDcorUP.setLegend(new DataTag[]{dsCorMSD.getTag()}, "MSD");
        sim.integrator.getEventManager().addListener(pumpMSDcor);
        DataPumpListener pumpPcor = new DataPumpListener(dsCorP, plotMSDcorUP.getDataSet().makeDataSink(), 1000);
        plotMSDcorUP.setLegend(new DataTag[]{dsCorP.getTag()}, "P");
        sim.integrator.getEventManager().addListener(pumpPcor);
        if (sim.potentialChoice != SimGlass.PotentialChoice.HS) {
            DataPumpListener pumpKEcor = new DataPumpListener(dsCorKE, plotMSDcorUP.getDataSet().makeDataSink(), 1000);
            plotMSDcorUP.setLegend(new DataTag[]{dsCorKE.getTag()}, "KE");
            sim.integrator.getEventManager().addListener(pumpKEcor);
        }

        if (dsMSDcorU != null) {
            meterMSD.addMSDSink(dsMSDcorU);
            plotMSDcorUP.setLegend(new DataTag[]{dsMSDcorU.getTag()}, "U");
            DataPumpListener pumpMSDcorU = new DataPumpListener(dsMSDcorU, plotMSDcorUP.getDataSet().makeDataSink(), 1000);
            sim.integrator.getEventManager().addListener(pumpMSDcorU);
        }

        meterMSD.addMSDSink(dsPMSDhistory);
        DisplayPlotXChart plotPMSDscatter = new DisplayPlotXChart();
        plotPMSDscatter.setLabel("MSD vs. P");
        plotPMSDscatter.setXLabel("P");
        DataPumpListener pumpPMSDhistory = new DataPumpListener(dsPMSDhistory, plotPMSDscatter.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpPMSDhistory);
        plotPMSDscatter.setDoDrawLines(new DataTag[]{dsPMSDhistory.getTag()}, false);
        plotPMSDscatter.getPlot().setYLabel("MSD");
        plotPMSDscatter.setDoLegend(false);


        DisplayPlotXChart plotHistogramP = new DisplayPlotXChart();
        DataPumpListener pumpHistogramP = new DataPumpListener(dsHistogramP, plotHistogramP.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpHistogramP);
        plotHistogramP.setXLabel("P");
        plotHistogramP.setDoLegend(false);
        plotHistogramP.getPlot().setYLog(true);

        DataSourceHistogramMSD dsHistogramMSD = new DataSourceHistogramMSD(sim.timeSource);
        meterMSD.addMSDSink(dsHistogramMSD);
        DisplayPlotXChart plotHistogramMSD = new DisplayPlotXChart();
        DataPumpListener pumpHistogramMSD = new DataPumpListener(dsHistogramMSD, plotHistogramMSD.getDataSet().makeDataSink(), 1000);
        sim.integrator.getEventManager().addListener(pumpHistogramMSD);
        plotHistogramMSD.setXLabel("MSD");
        plotHistogramMSD.setDoLegend(false);
        plotHistogramMSD.getPlot().setYLog(true);

        List<MeterStructureFactor.AtomSignalSourceByType> signalByTypes = new ArrayList<>();
        DeviceButtonGroup sfacButtons = null;
        int n = sim.box.getLeafList().size();
        double cut1 = 10;
        if (n > 500) cut1 /= Math.pow(n / 500.0, 1.0 / sim.getSpace().D());
        MeterStructureFactor meterSFac = new MeterStructureFactor(sim.box, cut1);
        meterSFac.setNormalizeByN(true);
        DataDump dumpSFac = new DataDump();
        DataFork forkSFac = new DataFork();
        signalByTypes.add((MeterStructureFactor.AtomSignalSourceByType) meterSFac.getSignalSource());
        DataPumpListener pumpSFac = new DataPumpListener(meterSFac, forkSFac, 1000);
        AccumulatorAverageFixed accSFac = new AccumulatorAverageFixed(1);  // just average, no uncertainty
        accSFac.setPushInterval(1);
        forkSFac.addDataSink(accSFac);
        forkSFac.addDataSink(dumpSFac);
        sim.integrator.getEventManager().addListener(pumpSFac);
        dataStreamPumps.add(pumpSFac);

        DisplayPlotXChart plotSFac = new DisplayPlotXChart();
        accSFac.addDataSink(plotSFac.makeSink("sfac"), new AccumulatorAverage.StatType[]{accSFac.AVERAGE});
        plotSFac.setLabel("SFac");
        plotSFac.getSeries("sfac")
                .setXYSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter)
                .setLabel("density");
        plotSFac.setYLog(true);
        plotSFac.getDataSet().setUpdatingOnAnyChange(true);

        double L = sim.box.getBoundary().getBoxSize().getX(0);

        MeterStructureFactor[] meterSFacMobility = new MeterStructureFactor[30];
        DataDump[] dumpSFacMobility = new DataDump[30];
        AccumulatorAverageFixed[] accSFacMobility = new AccumulatorAverageFixed[30];
        for (int i = 0; i < 30; i++) {
            AtomSignalMobility signalMobility = new AtomSignalMobility(configStorage);
            signalMobility.setPrevConfig(i + 1);
            meterSFacMobility[i] = new MeterStructureFactor(sim.box, 3, signalMobility);
            meterSFacMobility[i].setNormalizeByN(true);
            DataFork forkSFacMobility = new DataFork();
            accSFacMobility[i] = new AccumulatorAverageFixed(1);  // just average, no uncertainty
            accSFacMobility[i].setPushInterval(1);
            DataPump pumpSFacMobility = new DataPump(meterSFacMobility[i], forkSFacMobility);
            forkSFacMobility.addDataSink(accSFacMobility[i]);
            ConfigurationStoragePumper cspMobility = new ConfigurationStoragePumper(pumpSFacMobility, configStorage);
            cspMobility.setPrevStep(Math.max(i + 1, 7));
            configStorage.addListener(cspMobility);
            dataStreamPumps.add(pumpSFacMobility);
            dumpSFacMobility[i] = new DataDump();
            accSFacMobility[i].addDataSink(dumpSFacMobility[i], new AccumulatorAverage.StatType[]{accSFacMobility[i].AVERAGE});
        }

        MeterStructureFactor[] meterSFacMotion = new MeterStructureFactor[30];
        DataDump[] dumpSFacMotion = new DataDump[30];
        AccumulatorAverageFixed[] accSFacMotion = new AccumulatorAverageFixed[30];
        for (int i = 0; i < 30; i++) {
            AtomSignalMotion signalMotion = new AtomSignalMotion(configStorage, 0);
            signalMotion.setPrevConfig(i + 1);
            meterSFacMotion[i] = new MeterStructureFactor(sim.box, 3, signalMotion);
            meterSFacMotion[i].setNormalizeByN(true);
            DataFork forkSFacMotion = new DataFork();
            accSFacMotion[i] = new AccumulatorAverageFixed(1);  // just average, no uncertainty
            accSFacMotion[i].setPushInterval(1);
            DataPump pumpSFacMotion = new DataPump(meterSFacMotion[i], forkSFacMotion);
            forkSFacMotion.addDataSink(accSFacMotion[i]);
            ConfigurationStoragePumper cspMotion = new ConfigurationStoragePumper(pumpSFacMotion, configStorage);
            cspMotion.setPrevStep(Math.max(i + 1, 7));
            configStorage.addListener(cspMotion);
            dataStreamPumps.add(pumpSFacMotion);
            dumpSFacMotion[i] = new DataDump();
            accSFacMotion[i].addDataSink(dumpSFacMotion[i], new AccumulatorAverage.StatType[]{accSFacMotion[i].AVERAGE});
        }


        MeterFromDumps sfacFromDumps = new MeterFromDumps(dumpSFacMobility);
        DataPumpListener pumpSfacMobility = new DataPumpListener(sfacFromDumps, plotSFac.makeSink("mobility"), 500);
        sim.integrator.getEventManager().addListener(pumpSfacMobility);
        plotSFac.getSeries("mobility")
                .setXYSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter)
                .setLabel("mobility");

        MeterFromDumps sfacFromDumpsMotion = new MeterFromDumps(dumpSFacMotion);
        DataPumpListener pumpSfacMotion = new DataPumpListener(sfacFromDumpsMotion, plotSFac.makeSink("xmotion"), 500);
        sim.integrator.getEventManager().addListener(pumpSfacMotion);
        plotSFac.getSeries("xmotion")
                .setXYSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter)
                .setLabel("x motion");

        sfacButtons = new DeviceButtonGroup(sim.getController(), 5);
        sfacButtons.setLabel("B signal");
        AtomType typeB = sim.speciesB.getLeafType();
        sfacButtons.addButton("+1", new SFacButtonAction(signalByTypes, accSFac, typeB, +1));
        sfacButtons.addButton("-1", new SFacButtonAction(signalByTypes, accSFac, typeB, -1));
        double vB = sim.getSpace().powerD(sim.sigmaB);
        sfacButtons.addButton("+v", new SFacButtonAction(signalByTypes, accSFac, typeB, vB));
        sfacButtons.addButton("-v", new SFacButtonAction(signalByTypes, accSFac, typeB, -vB));
        sfacButtons.addButton("0", new SFacButtonAction(signalByTypes, accSFac, typeB, 0));
        sfacButtons.setSelected("+v");

        Vector[] wv = meterSFac.getWaveVectors();
        List<Vector> myWV = new ArrayList<>();
        double wvMax2 = 2.01 * Math.PI / L;
        for (Vector vector : wv) {
            int nd = 0;
            for (int i = 0; i < vector.getD(); i++) if (vector.getX(i) != 0) nd++;
            if (vector.squared() > wvMax2 * wvMax2 || nd > 1) continue;
            myWV.add(vector);
        }
        int minIntervalSfac2 = 8;
        wv = myWV.toArray(new Vector[0]);
        MeterStructureFactor[] meterSFacMotion2 = new MeterStructureFactor[30];
        int[] motionMap = StructureFactorComponentCorrelation.makeWaveVectorMap(wv, 0);
        StructureFactorComponentCorrelation sfcMotionCor = new StructureFactorComponentCorrelation(motionMap, configStorage);
        sfcMotionCor.setMinInterval(minIntervalSfac2);
        MeterStructureFactor[] meterSFacMobility2 = new MeterStructureFactor[30];
        int[] mobilityMap = StructureFactorComponentCorrelation.makeWaveVectorMap(wv, -1);
        StructureFactorComponentCorrelation sfcMobilityCor = new StructureFactorComponentCorrelation(mobilityMap, configStorage);
        sfcMobilityCor.setMinInterval(minIntervalSfac2);

        MeterStructureFactor meterSFacDensity2 = new MeterStructureFactor(sim.box, 3);
        meterSFacDensity2.setNormalizeByN(true);
        meterSFacDensity2.setWaveVec(wv);
        StructureFactorComponentCorrelation sfcDensityCor = new StructureFactorComponentCorrelation(mobilityMap, configStorage);
        sfcDensityCor.setMinInterval(0);
        DataSinkBlockAveragerSFac dsbaSfacDensity2 = new DataSinkBlockAveragerSFac(configStorage, 0, meterSFacDensity2);
        dsbaSfacDensity2.addSink(sfcDensityCor);
        DataPump pumpSFacDensity2 = new DataPump(meterSFacDensity2, dsbaSfacDensity2);
        ConfigurationStoragePumper cspDensity2 = new ConfigurationStoragePumper(pumpSFacDensity2, configStorage);
        configStorage.addListener(cspDensity2);
        cspDensity2.setPrevStep(0);
        DataSourceCorrelation dsCorSFacDensityMobility = new DataSourceCorrelation(configStorage, mobilityMap.length);
        dsbaSfacDensity2.addSink(dsCorSFacDensityMobility.makeReceiver(0));

        AtomSignalKineticEnergy atomSignalKE = new AtomSignalKineticEnergy();
        atomSignalKE.setDoSubtractAvg(1.5 * sim.integrator.getTemperature());
        MeterStructureFactor meterSFacKE = new MeterStructureFactor(sim.box, 3, atomSignalKE);
        meterSFacKE.setNormalizeByN(true);
        meterSFacKE.setWaveVec(wv);
        StructureFactorComponentCorrelation sfcKECor = new StructureFactorComponentCorrelation(mobilityMap, configStorage);
        sfcKECor.setMinInterval(0);
        DataSinkBlockAveragerSFac dsbaSfacKE = new DataSinkBlockAveragerSFac(configStorage, 0, meterSFacKE);
        dsbaSfacKE.addSink(sfcKECor);
        DataPump pumpSFacKECor = new DataPump(meterSFacKE, dsbaSfacKE);
        ConfigurationStoragePumper cspKE = new ConfigurationStoragePumper(pumpSFacKECor, configStorage);
        configStorage.addListener(cspKE);
        cspKE.setPrevStep(0);
        DataSourceCorrelation dsCorSFacKEMobility = new DataSourceCorrelation(configStorage, mobilityMap.length);
        dsbaSfacKE.addSink(dsCorSFacKEMobility.makeReceiver(0));

        AtomStressSource stressSource = null;
        if (sim.potentialChoice == SimGlass.PotentialChoice.HS) {
            AtomHardStressCollector ahsc = new AtomHardStressCollector(sim.timeSource);
            ((IntegratorHard) sim.integrator).addCollisionListener(ahsc);
            stressSource = ahsc;
        } else {
            PotentialCalculationForceSumGlass pcForce = new PotentialCalculationForceSumGlass(sim.box);
            ((IntegratorVelocityVerlet) sim.integrator).setForceSum(pcForce);
            stressSource = pcForce;
        }
        int[][] normalComps = new int[sim.getSpace().getD()][2];
        for (int i=0; i<normalComps.length; i++) {
            normalComps[i][0] = normalComps[i][1] = i;
        }
        AtomSignalStress signalStressNormal = new AtomSignalStress(stressSource, normalComps);

        MeterStructureFactor meterSFacStress2 = new MeterStructureFactor(sim.box, 3, signalStressNormal);
        meterSFacStress2.setNormalizeByN(true);
        meterSFacStress2.setWaveVec(wv);
        StructureFactorComponentCorrelation sfcStress2Cor = new StructureFactorComponentCorrelation(mobilityMap, configStorage);
        sfcStress2Cor.setMinInterval(0);
        DataSinkBlockAveragerSFac dsbaSfacStress2 = new DataSinkBlockAveragerSFac(configStorage, 0, meterSFacStress2);
        dsbaSfacStress2.addSink(sfcStress2Cor);
        DataPump pumpSFacStress2Cor = new DataPump(meterSFacStress2, dsbaSfacStress2);
        ConfigurationStoragePumper cspStress2 = new ConfigurationStoragePumper(pumpSFacStress2Cor, configStorage);
        configStorage.addListener(cspStress2);
        cspStress2.setPrevStep(0);
        DataSourceCorrelation dsCorSFacStress2Mobility = new DataSourceCorrelation(configStorage, mobilityMap.length);
        dsbaSfacStress2.addSink(dsCorSFacStress2Mobility.makeReceiver(0));

        DataSourceCorrelation dsCorSFacPackingMobility = new DataSourceCorrelation(configStorage, mobilityMap.length);
        DataSourceCorrelation dsCorSFacPackingDensity = new DataSourceCorrelation(configStorage, mobilityMap.length);
        StructureFactorComponentCorrelation sfcPackingCor = new StructureFactorComponentCorrelation(mobilityMap, configStorage);
        DataSourceCorrelation dsCorSFacDensityAMobility = new DataSourceCorrelation(configStorage, mobilityMap.length);
        DataSourceCorrelation dsCorSFacDensityADensity = new DataSourceCorrelation(configStorage, mobilityMap.length);
        StructureFactorComponentCorrelation sfcDensityACor = new StructureFactorComponentCorrelation(mobilityMap, configStorage);
        if (sim.potentialChoice == SimGlass.PotentialChoice.HS) {
            MeterStructureFactor.AtomSignalSourceByType atomSignalPacking = new MeterStructureFactor.AtomSignalSourceByType();
            atomSignalPacking.setAtomTypeFactor(sim.speciesB.getLeafType(), vB);
            MeterStructureFactor meterSFacPacking2 = new MeterStructureFactor(sim.box, 3, atomSignalPacking);
            meterSFacPacking2.setNormalizeByN(true);
            meterSFacPacking2.setWaveVec(wv);
            sfcPackingCor.setMinInterval(0);
            DataSinkBlockAveragerSFac dsbaSfacPacking2 = new DataSinkBlockAveragerSFac(configStorage, 0, meterSFacPacking2);
            dsbaSfacPacking2.addSink(sfcPackingCor);
            DataPump pumpSFacPacking2 = new DataPump(meterSFacPacking2, dsbaSfacPacking2);
            ConfigurationStoragePumper cspPacking2 = new ConfigurationStoragePumper(pumpSFacPacking2, configStorage);
            configStorage.addListener(cspPacking2);
            cspPacking2.setPrevStep(0);
            dsbaSfacPacking2.addSink(dsCorSFacPackingMobility.makeReceiver(0));
            dsbaSfacPacking2.addSink(dsCorSFacPackingDensity.makeReceiver(0));
            dsbaSfacDensity2.addSink(dsCorSFacPackingDensity.makeReceiver(1));
        } else {
            MeterStructureFactor.AtomSignalSourceByType atomSignalDensityA = new MeterStructureFactor.AtomSignalSourceByType();
            atomSignalDensityA.setAtomTypeFactor(sim.speciesB.getLeafType(), 0);
            MeterStructureFactor meterSFacDensityA2 = new MeterStructureFactor(sim.box, 3, atomSignalDensityA);
            meterSFacDensityA2.setNormalizeByN(true);
            meterSFacDensityA2.setWaveVec(wv);
            sfcDensityACor.setMinInterval(0);
            DataSinkBlockAveragerSFac dsbaSfacDensityA2 = new DataSinkBlockAveragerSFac(configStorage, 0, meterSFacDensityA2);
            dsbaSfacDensityA2.addSink(sfcDensityACor);
            DataPump pumpSFacDensityA2 = new DataPump(meterSFacDensityA2, dsbaSfacDensityA2);
            ConfigurationStoragePumper cspDensityA2 = new ConfigurationStoragePumper(pumpSFacDensityA2, configStorage);
            configStorage.addListener(cspDensityA2);
            cspDensityA2.setPrevStep(0);
            dsbaSfacDensityA2.addSink(dsCorSFacDensityAMobility.makeReceiver(0));
            dsbaSfacDensityA2.addSink(dsCorSFacDensityADensity.makeReceiver(0));
            dsbaSfacDensity2.addSink(dsCorSFacDensityADensity.makeReceiver(1));
        }

        for (int i = 0; i < 30; i++) {
            AtomSignalMotion signalMotion = new AtomSignalMotion(configStorage, 0);
            signalMotion.setPrevConfig(i + 1);
            meterSFacMotion2[i] = new MeterStructureFactor(sim.box, 3, signalMotion);
            meterSFacMotion2[i].setNormalizeByN(true);
            meterSFacMotion2[i].setWaveVec(wv);
            DataPump pumpSFacMotion2 = new DataPump(meterSFacMotion2[i], sfcMotionCor.makeSink(i, meterSFacMotion2[i]));
            ConfigurationStoragePumper cspMotion2 = new ConfigurationStoragePumper(pumpSFacMotion2, configStorage);
            cspMotion2.setPrevStep(i);
            cspMotion2.setBigStep(minIntervalSfac2);
            configStorage.addListener(cspMotion2);

            AtomSignalMobility signalMobility = new AtomSignalMobility(configStorage);
            signalMobility.setPrevConfig(i + 1);
            meterSFacMobility2[i] = new MeterStructureFactor(sim.box, 3, signalMobility);
            meterSFacMobility2[i].setNormalizeByN(true);
            meterSFacMobility2[i].setWaveVec(wv);
            DataFork sfacMobility2Fork = new DataFork();
            sfacMobility2Fork.addDataSink(sfcMobilityCor.makeSink(i, meterSFacMobility2[i]));
            DataPump pumpSFacMobility2 = new DataPump(meterSFacMobility2[i], sfacMobility2Fork);
            ConfigurationStoragePumper cspMobility2 = new ConfigurationStoragePumper(pumpSFacMobility2, configStorage);
            cspMobility2.setPrevStep(i);
            cspMobility2.setBigStep(minIntervalSfac2);
            configStorage.addListener(cspMobility2);
            sfacMobility2Fork.addDataSink(new StructorFactorComponentExtractor(meterSFacMobility2, i, dsCorSFacDensityMobility));
            if (sim.potentialChoice == SimGlass.PotentialChoice.HS) {
                sfacMobility2Fork.addDataSink(new StructorFactorComponentExtractor(meterSFacMobility2, i, dsCorSFacPackingMobility));
            } else {
                sfacMobility2Fork.addDataSink(new StructorFactorComponentExtractor(meterSFacMobility2, i, dsCorSFacDensityAMobility));
            }
            sfacMobility2Fork.addDataSink(new StructorFactorComponentExtractor(meterSFacMobility2, i, dsCorSFacKEMobility));
            sfacMobility2Fork.addDataSink(new StructorFactorComponentExtractor(meterSFacMobility2, i, dsCorSFacStress2Mobility));
        }

        DisplayPlotXChart plotSFacCor = new DisplayPlotXChart();
        plotSFacCor.getPlot().setTitle("correlation");
        plotSFacCor.setXLog(true);

        double fac = L / (2 * Math.PI);
        int[] foo = new int[mobilityMap.length];
        for (int j = 0; j < mobilityMap.length; j++) {
            if (foo[mobilityMap[j]] != 0) continue;
            foo[mobilityMap[j]] = 1;
            String label = String.format("%d,%d,%d", Math.round(Math.abs(myWV.get(j).getX(0)) * fac),
                    Math.round(Math.abs(myWV.get(j).getX(1)) * fac),
                    Math.round(Math.abs(myWV.get(j).getX(2)) * fac));

            StructureFactorComponentCorrelation.Meter m = sfcMobilityCor.makeMeter(mobilityMap[j]);
            DataPumpListener pumpSFacMobilityCor = new DataPumpListener(m, plotSFacCor.makeSink("mobility" + label), 1000);
            sim.integrator.getEventManager().addListener(pumpSFacMobilityCor);
            plotSFacCor.getSeries("mobility"+label)
                    .setLabel("mobility " + label);

            m = sfcDensityCor.makeMeter(mobilityMap[j]);
            DataPumpListener pumpSFacDensityCor = new DataPumpListener(m, plotSFacCor.makeSink("density"+label), 1000);
            sim.integrator.getEventManager().addListener(pumpSFacDensityCor);
            plotSFacCor.getSeries("density"+label)
                    .setLabel("density " + label);

            m = sfcKECor.makeMeter(mobilityMap[j]);
            DataPumpListener pumpSFacKEsfcCor = new DataPumpListener(m, plotSFacCor.makeSink("KE" + label), 1000);
            sim.integrator.getEventManager().addListener(pumpSFacKEsfcCor);

            plotSFacCor.getSeries("KE"+label).setLabel("KE " + label);

            m = sfcStress2Cor.makeMeter(mobilityMap[j]);
            DataPumpListener pumpSFacStress2sfcCor = new DataPumpListener(m, plotSFacCor.makeSink("stress" + label), 1000);
            sim.integrator.getEventManager().addListener(pumpSFacStress2sfcCor);
            plotSFacCor.getSeries("stress"+label).setLabel("stress " + label);

            if (sim.potentialChoice == SimGlass.PotentialChoice.HS) {
                m = sfcPackingCor.makeMeter(mobilityMap[j]);
                DataPumpListener pumpSFacPackingCor = new DataPumpListener(m, plotSFacCor.makeSink("packing"+label), 1000);
                sim.integrator.getEventManager().addListener(pumpSFacPackingCor);
                plotSFacCor.getSeries("packing"+label).setLabel("packing " + label);
            } else {
                m = sfcDensityACor.makeMeter(mobilityMap[j]);
                DataPumpListener pumpSFacDensityACor = new DataPumpListener(m, plotSFacCor.makeSink("densityA"+label), 1000);
                sim.integrator.getEventManager().addListener(pumpSFacDensityACor);
                plotSFacCor.getSeries("densityA"+label).setLabel("densityA " + label);
            }
            DataSourceCorrelation.Meter mm = dsCorSFacDensityMobility.makeMeter(j);
            DataPumpListener pumpCorSFacDensityMobility = new DataPumpListener(mm, plotSFacCor.makeSink("d-m"+label), 1000);
            sim.integrator.getEventManager().addListener(pumpCorSFacDensityMobility);
            plotSFacCor.getSeries("d-m"+label).setLabel("d-m " + label);

            if (sim.potentialChoice == SimGlass.PotentialChoice.HS) {
                mm = dsCorSFacPackingMobility.makeMeter(j);
                DataPumpListener pumpCorSFacPackingMobility = new DataPumpListener(mm, plotSFacCor.makeSink("p-m"+label), 1000);
                sim.integrator.getEventManager().addListener(pumpCorSFacPackingMobility);
                plotSFacCor.getSeries("p-m"+label).setLabel("p-m " + label);

                mm = dsCorSFacPackingDensity.makeMeter(j);
                DataPumpListener pumpCorSFacPackingDensity = new DataPumpListener(mm, plotSFacCor.makeSink("p-d"+label), 1000);
                sim.integrator.getEventManager().addListener(pumpCorSFacPackingDensity);
                plotSFacCor.getSeries("p-d"+label).setLabel("p-d " + label);
            } else {
                mm = dsCorSFacDensityAMobility.makeMeter(j);
                DataPumpListener pumpCorSFacDensityAMobility = new DataPumpListener(mm, plotSFacCor.makeSink("da-m"+label), 1000);
                sim.integrator.getEventManager().addListener(pumpCorSFacDensityAMobility);
                plotSFacCor.getSeries("da-m"+label).setLabel("da-m " + label);

                mm = dsCorSFacDensityADensity.makeMeter(j);
                DataPumpListener pumpCorSFacDensityADensity = new DataPumpListener(mm, plotSFacCor.makeSink("da-d"+label), 1000);
                sim.integrator.getEventManager().addListener(pumpCorSFacDensityADensity);
                plotSFacCor.getSeries("da-d"+label).setLabel("da-d " + label);
            }
            mm = dsCorSFacKEMobility.makeMeter(j);
            DataPumpListener pumpCorSFacKEMobility = new DataPumpListener(mm, plotSFacCor.makeSink("KE-m"+label), 1000);
            sim.integrator.getEventManager().addListener(pumpCorSFacKEMobility);
            plotSFacCor.getSeries("KE-m"+label).setLabel("KE-m " + label);

            mm = dsCorSFacStress2Mobility.makeMeter(j);
            DataPumpListener pumpCorSFacStress2Mobility = new DataPumpListener(mm, plotSFacCor.makeSink("s-m"+label), 1000);
            sim.integrator.getEventManager().addListener(pumpCorSFacStress2Mobility);
            plotSFacCor.getSeries("s-m"+label).setLabel("s-m " + label);
        }

        foo = new int[motionMap.length];
        for (int j = 0; j < motionMap.length; j++) {
            if (foo[motionMap[j]] != 0) continue;
            foo[motionMap[j]] = 1;
            String label = String.format("%d,%d,%d", Math.round(myWV.get(j).getX(0) * fac),
                    Math.round(myWV.get(j).getX(1) * fac),
                    Math.round(myWV.get(j).getX(2) * fac));

            StructureFactorComponentCorrelation.Meter m = sfcMotionCor.makeMeter(motionMap[j]);
            DataPumpListener pumpSFacMotionCor = new DataPumpListener(m, plotSFacCor.makeSink("motion"+label), 1000);
            sim.integrator.getEventManager().addListener(pumpSFacMotionCor);
            plotSFacCor.getSeries("motion"+label).setLabel("motion " + label);
        }

        DeviceSlider sfacPrevConfig = new DeviceSlider(sim.getController(), new Modifier() {
            @Override
            public void setValue(double newValue) {
                int idx = (int) Math.round(newValue);
                sfacFromDumps.setDumpIndex(idx);
                pumpSfacMobility.actionPerformed();
                sfacFromDumpsMotion.setDumpIndex(idx);
                pumpSfacMotion.actionPerformed();
            }

            @Override
            public double getValue() {
                return sfacFromDumps.getDumpIndex();
            }

            @Override
            public Dimension getDimension() {
                return Null.DIMENSION;
            }

            @Override
            public String getLabel() {
                return "prev config";
            }
        });
        sfacPrevConfig.setMinimum(0);
        sfacPrevConfig.setMaximum(30);
        sfacPrevConfig.setNMajor(5);
        sfacPrevConfig.setValue(5);
        sfacPrevConfig.setLabel("Previous config (mobility)");
        sfacPrevConfig.setShowBorder(true);

        MeterStructureFactor meterSFacNormal = new MeterStructureFactor(sim.box, 3, signalStressNormal);
        AccumulatorAverageFixed accSFacNormalStress = new AccumulatorAverageFixed(1);
        DataPumpListener pumpSFacNormal = new DataPumpListener(meterSFacNormal, accSFacNormalStress, 100);
        accSFacNormalStress.setPushInterval(5);
        sim.integrator.getEventManager().addListener(pumpSFacNormal);
        dataStreamPumps.add(pumpSFacNormal);

        accSFacNormalStress.addDataSink(plotSFac.makeSink("normalstress"), new AccumulatorAverage.StatType[]{accSFac.AVERAGE});
        plotSFac.getSeries("normalstress")
                .setXYSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter)
                .setLabel("normal stress");

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
                    meterMSD.reset();
                    meterMSDA.reset();
                    meterMSDB.reset();
                    meterPerc.zeroData();
                    meterQ4.zeroData();
                    dsMSDcorP.setEnabled(false);
                    dsCorP.setEnabled(false);
                    if (dsCorKE!=null) dsCorKE.setEnabled(false);
                    dsCorMSD.setEnabled(false);
                    dsHistogramP.setEnabled(false);
                    dsPMSDhistory.setEnabled(false);
                    dpAutocor.reset();
                    dpxyAutocor.reset();
                    pTensorAccum.setEnabled(false);
                    pTensorAccum.reset();
                    dsMSDcorU.setEnabled(false);
                    meterFs.reset();
                    meterFsA.reset();
                    meterFsB.reset();
                    meterF.reset();
                    meterVAC.reset();
                    meterCorrelationAA.zeroData();
                    meterCorrelationAB.zeroData();
                    meterCorrelationBB.zeroData();
                    meterCorrelation.zeroData();
                    meterCorrelationPerp.zeroData();
                    meterCorrelationPar.zeroData();
                    meterCorrelationMag.zeroData();
                    meterCorrelationAAMag.zeroData();
                    meterCorrelationABMag.zeroData();
                    meterCorrelationBBMag.zeroData();
                    diameterHash.setFac(1.0);
                    accSFac.reset();
                    sfcMotionCor.reset();
                    sfcDensityACor.reset();
                    sfcDensityCor.reset();
                    sfcKECor.reset();
                    sfcMobilityCor.reset();
                    sfcPackingCor.reset();
                    dsCorSFacDensityMobility.reset();
                    dsCorSFacDensityADensity.reset();
                    dsCorSFacDensityAMobility.reset();
                    dsCorSFacKEMobility.reset();
                    dsCorSFacPackingDensity.reset();
                    dsCorSFacPackingMobility.reset();
                    for (AccumulatorAverageFixed accumulatorAverageFixed : accSFacMobility) {
                        accumulatorAverageFixed.reset();
                    }
                    for (AccumulatorAverageFixed accumulatorAverageFixed : accSFacMotion) {
                        accumulatorAverageFixed.reset();
                    }
                    accSFacNormalStress.reset();
                    accPerc0.reset();
                    meterCorrelationSelf.reset();
                    meterCorrelationSelfMagA.reset();
                    meterCorrelationSelfMagB.reset();
                } else {
                    dbox.setAtomTestDoDisplay(atomFilterDeviation);
                    sim.integrator.setIntegratorMC(null, 0);
                    sim.integrator.setIsothermal(false);
                    configStorageLinear.reset();
                    configStorageLinear.setEnabled(true);
                    configStorage.reset();
                    configStorage.setEnabled(true);
                    meterMSD.reset();
                    meterMSDA.reset();
                    meterMSDB.reset();
                    meterPerc.zeroData();
                    meterQ4.zeroData();
                    dsMSDcorP.setEnabled(true);
                    dsCorP.setEnabled(true);
                    dsCorKE.setEnabled(true);
                    dsCorMSD.setEnabled(true);
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
                    meterVAC.reset();
                    meterCorrelationAA.zeroData();
                    meterCorrelationAB.zeroData();
                    meterCorrelationBB.zeroData();
                    meterCorrelation.zeroData();
                    meterCorrelationPerp.zeroData();
                    meterCorrelationPar.zeroData();
                    meterCorrelationMag.zeroData();
                    meterCorrelationAAMag.zeroData();
                    meterCorrelationABMag.zeroData();
                    meterCorrelationBBMag.zeroData();
                    diameterHash.setFac(showDispCheckbox.getState() ? 0.5 : 1.0);
                    accSFac.reset();
                    sfcMotionCor.reset();
                    sfcDensityACor.reset();
                    sfcDensityCor.reset();
                    sfcKECor.reset();
                    sfcMobilityCor.reset();
                    sfcPackingCor.reset();
                    dsCorSFacDensityMobility.reset();
                    dsCorSFacDensityADensity.reset();
                    dsCorSFacDensityAMobility.reset();
                    dsCorSFacKEMobility.reset();
                    dsCorSFacPackingDensity.reset();
                    dsCorSFacPackingMobility.reset();
                    for (AccumulatorAverageFixed accumulatorAverageFixed : accSFacMobility) {
                        accumulatorAverageFixed.reset();
                    }
                    for (AccumulatorAverageFixed accumulatorAverageFixed : accSFacMotion) {
                        accumulatorAverageFixed.reset();
                    }
                    accSFacNormalStress.reset();
                    accPerc0.reset();
                    meterCorrelationSelf.reset();
                    meterCorrelationSelfMagA.reset();
                    meterCorrelationSelfMagB.reset();
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
        getPanel().tabbedPane.add("P cor MSD", plotsPane);

        JPanel plotSFacPanel = new JPanel(new MigLayout("fill"));
        plotSFacPanel.add(sfacButtons.graphic(), "alignx center");
        plotSFacPanel.add(sfacPrevConfig.graphic(), "wrap, alignx center");
        plotSFacPanel.add(plotSFac.graphic(), "grow");
        plotSFacPanel.add(plotSFacCor.graphic(), "grow");
        JPanel sfacWidgetPanel = new JPanel();

        plotSFacPanel.add(sfacWidgetPanel);

        JScrollPane plotsPaneSFac = new JScrollPane(plotSFacPanel);
        getPanel().tabbedPane.add("SFac", plotsPaneSFac);
    }

    public static void main(String[] args) {
        SimGlass.GlassParams params = new SimGlass.GlassParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.potential = SimGlass.PotentialChoice.HS;
            params.nA = 250;
            params.nB = 250;
            params.density = 1.52;
            params.D = 3;
        }
        SimGlass sim = new SimGlass(params.D, params.nA, params.nB, params.density, params.temperature, params.doSwap, params.potential, params.tStep);

        GlassGraphic ljmdGraphic = new GlassGraphic(sim);
        SimulationGraphic.makeAndDisplayFrame
                (ljmdGraphic.getPanel(), sim.potentialChoice + APP_NAME);
        if (params.potential == SimGlass.PotentialChoice.HS) {
            JFrame f = new JFrame();
            f.setSize(700, 500);
            JPanel panel = new JPanel();
            panel.add(new JLabel("<html><div style='width: 200px;'>Note: high-density simulations will be slightly delayed in starting up, as the simulation works to generate a configuration without overlaps. Controls probably won't work properly until this dialog disappears. </div></html>"));
            f.getContentPane().add(panel);
            f.pack();
            f.setTitle("Generating configuration");
            f.setLocationRelativeTo(null);
            f.setVisible(true);
            sim.getController().addActivity(sim.makeInitConfigActivity()).future.whenComplete((res, ex) -> {
                ljmdGraphic.getController().getResetAveragesButton().getAction().actionPerformed();
                f.dispatchEvent(new WindowEvent(f, WindowEvent.WINDOW_CLOSING));
            });
        }
        sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
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
        protected final java.util.List<MeterStructureFactor.AtomSignalSourceByType> signalSources;
        protected final AccumulatorAverageFixed acc;

        public SFacButtonAction(java.util.List<MeterStructureFactor.AtomSignalSourceByType> signalSources, AccumulatorAverageFixed acc, AtomType type, double value) {
            this.signalSources = signalSources;
            this.acc = acc;
            this.type = type;
            this.value = value;
        }

        public void actionPerformed() {
            for (MeterStructureFactor.AtomSignalSourceByType s : signalSources) {
                s.setAtomTypeFactor(type, value);
            }
            acc.reset();
        }
    }

    public static class MeterFromDumps implements IDataSource {
        private final DataDump[] dumps;
        private int dumpIndex;
        private final DataTag tag;
        private IDataInfo dataInfo;

        public MeterFromDumps(DataDump[] dumps) {
            this.dumps = dumps;
            tag = new DataTag();
            setDumpIndex(0);
        }

        public void setDumpIndex(int newDumpIndex) {
            dumpIndex = newDumpIndex;
            dataInfo = dumps[dumpIndex].getDataInfo().getFactory().makeDataInfo();
            dataInfo.addTag(tag);
        }

        public int getDumpIndex() {
            return dumpIndex;
        }

        public DataTag getTag() {
            return tag;
        }

        public IData getData() {
            return dumps[dumpIndex].getData();
        }

        public IDataInfo getDataInfo() {
            return dataInfo;
        }
    }
}
