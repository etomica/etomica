/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic.overlap;

import etomica.action.IAction;
import etomica.atom.DiameterHashByType;
import etomica.data.*;
import etomica.data.DataSplitter.IDataSinkFactory;
import etomica.data.histogram.HistogramExpanding;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.history.HistoryComplete;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.data.types.DataGroup;
import etomica.graphics.*;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.math.function.Function;
import etomica.math.function.IFunction;
import etomica.math.numerical.AkimaSpline;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.overlap.AlphaSource;
import etomica.overlap.DataOverlap;
import etomica.overlap.DataOverlap.DataSourceOverlapAvg;
import etomica.overlap.DataOverlap.DataSourceOverlapAvgCollapsing;
import etomica.overlap.DataOverlap.DataSourceOverlapLogAvg;
import etomica.overlap.IntegratorOverlap;
import etomica.overlap.MeterOverlap;
import etomica.space.Space;
import etomica.space1d.Vector1D;
import etomica.units.Pixel;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.*;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;

public class MultiharmonicGraphicMC extends SimulationGraphic {

    private final static String APP_NAME = "MultiharmonicMC";
    protected final static int REPAINT_INTERVAL = 300000;
    protected DataSplitter splitterA, splitterB;
    protected DataSourceOverlapAvg overlapAvgA, overlapAvgB;
    protected DataSourceAlphaFE dataSourceRefFE, dataSourceTargetFE;
    protected DataSourceAlphaFE dataSourceRefStdev, dataSourceTargetStdev;
    protected DataSourceAlphaFE dataSourceRefBias, dataSourceTargetBias;
    protected double alphaChoice = 1;
    protected final MultiharmonicMC sim;
    protected boolean makeDWPlot = false;

    /**
     * 
     */
    public MultiharmonicGraphicMC(MultiharmonicMC simulation, Space _space) {
        super(simulation, GRAPHIC_ONLY, APP_NAME, REPAINT_INTERVAL);
        this.sim = simulation;
        
        
        final MeterPotentialEnergy meterPEAinA = new MeterPotentialEnergy(sim.potentialMasterA, sim.boxA);
        final MeterPotentialEnergy meterPEAinB = new MeterPotentialEnergy(sim.potentialMasterA, sim.boxB);
        final MeterPotentialEnergy meterPEBinA = new MeterPotentialEnergy(sim.potentialMasterB, sim.boxA);
        final MeterPotentialEnergy meterPEBinB = new MeterPotentialEnergy(sim.potentialMasterB, sim.boxB);
        
        final int na = 11;
        final MeterOverlap meterOverlapA = new MeterOverlap(meterPEAinA, meterPEBinA, 1.0, true);
        meterOverlapA.setNumAlpha(na);
        meterOverlapA.setAlphaRange(1, 2*(na-1)*Math.log(2));
        final MeterOverlap meterOverlapB = new MeterOverlap(meterPEAinB, meterPEBinB, 1.0, false);
        meterOverlapB.setNumAlpha(na);
        meterOverlapB.setAlphaRange(1, 2*(na-1)*Math.log(2));

        final DataFork forkA = new DataFork();
        final DataFork forkB = new DataFork();
        
        
        DataPumpListener dataPumpA = new DataPumpListener(meterOverlapA, forkA);
        sim.integratorA.getEventManager().addListener(dataPumpA);

        DataPumpListener dataPumpB = new DataPumpListener(meterOverlapB, forkB);
        sim.integratorB.getEventManager().addListener(dataPumpB);

        final int dataInterval = 1;
        final boolean doCollapsing = true;
        final IDataSinkFactory dataSinkFactory = new IDataSinkFactory() {
            public IDataSink makeDataSink(int i) {
                AccumulatorAverageCollapsingLogAB acci = new AccumulatorAverageCollapsingLogAB();
                acci.setNumRawDataDoubles(8);
                return acci;
            }
        };
        if (doCollapsing) {
            
            splitterA = new DataSplitter();
            splitterA.setDataSinkFactory(dataSinkFactory);
            forkA.addDataSink(splitterA);
            for (int i=0; i<na; i++) {
                ((AccumulatorAverageCollapsingLogAB)splitterA.getDataSink(i)).setMaxSample(1);
            }
            overlapAvgA = new DataOverlap.DataSourceOverlapAvgCollapsingSplit(splitterA);
            
            splitterB = new DataSplitter();
            splitterB.setDataSinkFactory(dataSinkFactory);
            forkB.addDataSink(splitterB);
            for (int i=0; i<na; i++) {
                ((AccumulatorAverageCollapsingLogAB)splitterB.getDataSink(i)).setMaxSample(1.0/meterOverlapB.getAlpha(i));
            }
            overlapAvgB = new DataOverlap.DataSourceOverlapAvgCollapsingSplit(splitterB);
        }
        else {
            AccumulatorAverageFixed accFixedA = new AccumulatorAverageFixed();
            forkA.addDataSink(accFixedA);
            overlapAvgA = new DataOverlap.DataSourceOverlapAvgSimple(accFixedA);

            AccumulatorAverageFixed accFixedB = new AccumulatorAverageFixed();
            forkB.addDataSink(accFixedB);
            overlapAvgB = new DataOverlap.DataSourceOverlapAvgSimple(accFixedB);
        }
        DisplayPlot dWPlot = null;
        if (makeDWPlot) {
            dWPlot = new DisplayPlot();
            dWPlot.setLabel("dW");
            dWPlot.getPlot().setXLabel("Ui-UW");
            dWPlot.getDataSet().setUpdatingOnAnyChange(true);
            final AccumulatorHistogram dAW = new AccumulatorHistogram(new HistogramExpanding(0.05));
            forkA.addDataSink(new IDataSink() {
                final DataDouble myData = new DataDouble();
                final DataInfoDouble myDataInfo = new DataInfoDouble("dAW", Energy.DIMENSION);
                public void putDataInfo(IDataInfo inputDataInfo) {
                    dAW.putDataInfo(myDataInfo);
                }
                
                public void putData(IData inputData) {
                    myData.x = Math.log(inputData.getValue((inputData.getLength()-1)/2));
                    dAW.putData(myData);
                }
            });
            dAW.setPushInterval(100000);
            dAW.addDataSink(dWPlot.getDataSet().makeDataSink());
            dWPlot.setLegend(new DataTag[]{dAW.getTag()}, "dAW");
            final AccumulatorHistogram dBW = new AccumulatorHistogram(new HistogramExpanding(0.1));
            forkB.addDataSink(new IDataSink() {
                final DataDouble myData = new DataDouble();
                final DataInfoDouble myDataInfo = new DataInfoDouble("dBW", Energy.DIMENSION);
                public void putDataInfo(IDataInfo inputDataInfo) {
                    dBW.putDataInfo(myDataInfo);
                }
                
                public void putData(IData inputData) {
                    myData.x = Math.log(inputData.getValue((inputData.getLength()-1)/2));
                    dBW.putData(myData);
                }
            });
            dBW.setPushInterval(100000);
            dBW.addDataSink(dWPlot.getDataSet().makeDataSink());
            dWPlot.setLegend(new DataTag[]{dBW.getTag()}, "dBW");
        }
        
        final DataOverlap dsvo = new DataOverlapCaching(overlapAvgA, overlapAvgB, meterOverlapA, sim.integratorOS);
        sim.integratorOS.setReferenceFracSource(dsvo);
        sim.integratorOS.setAdjustStepFraction(true);
        sim.integratorOS.setNumSubSteps(1);
        sim.integratorOS.setAdjustInterval(dataInterval);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        sim.integratorOS.getEventManager().addListener(new IntegratorListener() {
            
            public void integratorStepStarted(IntegratorEvent e) {}
            public void integratorInitialized(IntegratorEvent e) {}
            
            public void integratorStepFinished(IntegratorEvent e) {
                long stepCount = sim.integratorOS.getStepCount();
                if (stepCount < 100) return;
                long oldInterval = sim.integratorOS.getAdjustInterval();
                if (stepCount % oldInterval == 0) {
                    int newInterval = (int)(stepCount / 50);
                    if (newInterval == oldInterval * 2) {
                        sim.integratorOS.setAdjustInterval((int)(stepCount / 50));
                    }
                }
            }
        });
        
        final DisplayBox displayBoxA = getDisplayBox(sim.boxA);
        remove(displayBoxA);
        final DisplayBox displayBoxB = getDisplayBox(sim.boxB);
        remove(displayBoxB);
        final Device controllerButtons = getController();
        getPanel().graphicsPanel.remove(controllerButtons.graphic());
        getPanel().footerPanel.add(controllerButtons.graphic());
        getPanel().graphicsPanel.setLayout(new GridBagLayout());
        
        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();
        dataStreamPumps.add(dataPumpA);

        displayBoxA.setPixelUnit(new Pixel(350/sim.boxA.getBoundary().getBoxSize().getX(0)));
        ((DiameterHashByType)displayBoxA.getDiameterHash()).setDiameter(sim.species.getLeafType(), 0.02);
        displayBoxB.setPixelUnit(new Pixel(350/sim.boxB.getBoundary().getBoxSize().getX(0)));
        displayBoxB.setDiameterHash(displayBoxA.getDiameterHash());

        final DisplayPlot fePlot = new DisplayPlot();
        final DisplayPlot feLnPlot = new DisplayPlot();
        feLnPlot.getPlot().setXLog(true);
        
        final DataSourceCountSteps stepCounter = new DataSourceCountSteps(sim.integratorOS);

        final ArrayList<DataPumpListener> allPumps = new ArrayList<DataPumpListener>();
        
        DataSourceScalar feAntibiasDataSource = new DataSourceScalar("antibias", Null.DIMENSION) {
            public double getDataAsScalar() {
                long refSteps = sim.integratorA.getStepCount();
                long targetSteps = sim.integratorB.getStepCount();
                double myAlpha = (dsvo.getOverlapAverageAndError()[0]*refSteps)/targetSteps;
                double refChi = dsvo.getAverageAndError(true, myAlpha)[0];
                double targetChi = dsvo.getAverageAndError(false, myAlpha)[0];
                long refCount = sim.integratorA.getStepCount();
                double antibiasRefChi = (refChi * refCount + 1) / (refCount+1);
                double lnRefDiff = 0.5*Math.log(refChi/antibiasRefChi);
                long targetCount = sim.integratorB.getStepCount();
                double antibiasTargetChi = (targetChi * targetCount + 1/myAlpha) / (targetCount+1);
                double lnTargetDiff = 0.5*(Math.log(targetChi/antibiasTargetChi));
                return -Math.log(refChi/targetChi) + lnRefDiff - lnTargetDiff; 
            }
        };
        
        DataSourceScalar feAntibiasPDataSource = new DataSourceScalar("antibias", Null.DIMENSION) {
            public double getDataAsScalar() {
                long refSteps = sim.integratorA.getStepCount();
                long targetSteps = sim.integratorB.getStepCount();
                double myAlpha = (dsvo.getOverlapAverageAndError()[0]*refSteps)/targetSteps;
                double myErr = dsvo.getLogAverageAndError(myAlpha)[1];
                double refChi = dsvo.getAverageAndError(true, myAlpha)[0];
                double targetChi = dsvo.getAverageAndError(false, myAlpha)[0];
                long refCount = sim.integratorA.getStepCount();
                double antibiasRefChi = (refChi * refCount + 1) / (refCount+1);
                double lnRefDiff = 0.5*Math.log(refChi/antibiasRefChi);
                long targetCount = sim.integratorB.getStepCount();
                double antibiasTargetChi = (targetChi * targetCount + 1/myAlpha) / (targetCount+1);
                double lnTargetDiff = 0.5*Math.log(targetChi/antibiasTargetChi);
                double lnDiffSum = Math.abs(lnRefDiff + lnTargetDiff);
                return -Math.log(refChi/targetChi) + lnRefDiff - lnTargetDiff + myErr + lnDiffSum; 
            }
        };

        DataSourceScalar feDataSource = new DataSourceScalar("average", Null.DIMENSION) {
            public double getDataAsScalar() {
                long refSteps = sim.integratorA.getStepCount();
                long targetSteps = sim.integratorB.getStepCount();
                double bestAlpha = (dsvo.getOverlapAverageAndError()[0]*refSteps)/targetSteps;
                return -dsvo.getLogAverageAndError(bestAlpha)[0];
            }
        };
        DataSourceScalar fePDataSource = new DataSourceScalar("error+", Null.DIMENSION) {
            public double getDataAsScalar() {
                long refSteps = sim.integratorA.getStepCount();
                long targetSteps = sim.integratorB.getStepCount();
                double bestAlpha = (dsvo.getOverlapAverageAndError()[0]*refSteps)/targetSteps;
                double[] avgErr = dsvo.getLogAverageAndError(bestAlpha);
                return -avgErr[0] + avgErr[1];
            }
        };

        AccumulatorHistory feLHistory = new AccumulatorHistory(new HistoryCollapsingDiscard(102, 3));
        DataPumpListenerSmart feLPump = new DataPumpListenerSmart(feDataSource, feLHistory, dataInterval, (HistoryCollapsingDiscard)feLHistory.getHistory());
        sim.integratorOS.getEventManager().addListener(feLPump);
        feLHistory.setTimeDataSource(stepCounter);
        feLHistory.setPushInterval(1);
        feLHistory.addDataSink(fePlot.getDataSet().makeDataSink());
        fePlot.setLegend(new DataTag[]{feLHistory.getTag()}, "fe");

        AccumulatorHistory feLHistoryLn = new AccumulatorHistory(new HistoryComplete());
        DataPumpListenerPowStrict feLPumpLn = new DataPumpListenerPowStrict(feDataSource, feLHistoryLn);
        sim.integratorOS.getEventManager().addListener(feLPumpLn);
        feLHistoryLn.setTimeDataSource(stepCounter);
        feLHistoryLn.setPushInterval(1);
        feLHistoryLn.addDataSink(feLnPlot.getDataSet().makeDataSink());
        feLnPlot.setLegend(new DataTag[]{feLHistoryLn.getTag()}, "fe");

        AccumulatorHistory feLPHistory = new AccumulatorHistory(new HistoryCollapsingDiscard(102, 3));
        DataPumpListenerSmart feLPPump = new DataPumpListenerSmart(fePDataSource, feLPHistory, dataInterval, (HistoryCollapsingDiscard)feLPHistory.getHistory());
        sim.integratorOS.getEventManager().addListener(feLPPump);
        feLPHistory.setTimeDataSource(stepCounter);
        feLPHistory.setPushInterval(1);
        feLPHistory.setDataSink(fePlot.getDataSet().makeDataSink());
        fePlot.setLegend(new DataTag[]{feLPHistory.getTag()}, "fe+");

        AccumulatorHistory feLPHistoryLn = new AccumulatorHistory(new HistoryComplete());
        DataPumpListenerPowStrict feLPPumpLn = new DataPumpListenerPowStrict(fePDataSource, feLPHistoryLn);
        sim.integratorOS.getEventManager().addListener(feLPPumpLn);
        feLPHistoryLn.setTimeDataSource(stepCounter);
        feLPHistoryLn.setPushInterval(1);
        feLPHistoryLn.setDataSink(feLnPlot.getDataSet().makeDataSink());
        feLnPlot.setLegend(new DataTag[]{feLPHistoryLn.getTag()}, "fe+");
        
        AccumulatorHistory feAntibiasHistory = new AccumulatorHistory(new HistoryCollapsingDiscard(102, 3));
        DataPumpListenerSmart feAntibiasPump = new DataPumpListenerSmart(feAntibiasDataSource, feAntibiasHistory, dataInterval, (HistoryCollapsingDiscard)feAntibiasHistory.getHistory());
        sim.integratorOS.getEventManager().addListener(feAntibiasPump);
        feAntibiasHistory.setTimeDataSource(stepCounter);
        feAntibiasHistory.addDataSink(fePlot.getDataSet().makeDataSink());
        fePlot.setLegend(new DataTag[]{feAntibiasHistory.getTag()}, "fe(ab)");
        
        AccumulatorHistory feAntibiasHistoryLn = new AccumulatorHistory(new HistoryComplete());
        DataPumpListenerPowStrict feAntibiasPumpLn = new DataPumpListenerPowStrict(feAntibiasDataSource, feAntibiasHistoryLn);
        sim.integratorOS.getEventManager().addListener(feAntibiasPumpLn);
        feAntibiasHistoryLn.setTimeDataSource(stepCounter);
        feAntibiasHistoryLn.addDataSink(feLnPlot.getDataSet().makeDataSink());
        feLnPlot.setLegend(new DataTag[]{feAntibiasHistoryLn.getTag()}, "fe(ab)");

        AccumulatorHistory feAntibiasErrHistory = new AccumulatorHistory(new HistoryCollapsingDiscard(102, 3));
        DataPumpListenerSmart feAntibiasErrPump = new DataPumpListenerSmart(feAntibiasPDataSource, feAntibiasErrHistory, dataInterval, (HistoryCollapsingDiscard)feAntibiasErrHistory.getHistory());
        sim.integratorOS.getEventManager().addListener(feAntibiasErrPump);
        feAntibiasErrHistory.setTimeDataSource(stepCounter);
        feAntibiasErrHistory.addDataSink(fePlot.getDataSet().makeDataSink());
        fePlot.setLegend(new DataTag[]{feAntibiasErrHistory.getTag()}, "fe(ab)+");

        AccumulatorHistory feAntibiasErrHistoryLn = new AccumulatorHistory(new HistoryComplete());
        DataPumpListenerPowStrict feAntibiasErrPumpLn = new DataPumpListenerPowStrict(feAntibiasPDataSource, feAntibiasErrHistoryLn);
        sim.integratorOS.getEventManager().addListener(feAntibiasErrPumpLn);
        feAntibiasErrHistoryLn.setTimeDataSource(stepCounter);
        feAntibiasErrHistoryLn.addDataSink(feLnPlot.getDataSet().makeDataSink());
        feLnPlot.setLegend(new DataTag[]{feAntibiasErrHistoryLn.getTag()}, "fe(ab)+");

        DataSourceScalar fracDataSource = new DataSourceScalar("frac A", Fraction.DIMENSION) {
            public double getDataAsScalar() {
                return sim.integratorOS.getRefStepFraction();
            }
        };
        AccumulatorHistory fracHistory = new AccumulatorHistory(new HistoryCollapsingDiscard(102, 3));
        fracHistory.setTimeDataSource(stepCounter);
        DataPumpListener fracPump = new DataPumpListener(fracDataSource, fracHistory, 1);
        sim.integratorOS.getEventManager().addListener(fracPump);
        DisplayPlot fracPlot = new DisplayPlot();
        fracHistory.setDataSink(fracPlot.getDataSet().makeDataSink());
        fracPlot.setLegend(new DataTag[]{fracHistory.getTag()}, "actual");
        DataSourceScalar frac2DataSource = new DataSourceScalar("frac A", Fraction.DIMENSION) {
            public double getDataAsScalar() {
                return dsvo.getIdealRefFraction(sim.integratorOS.getRefStepFraction());
            }
        };
        AccumulatorHistory frac2History = new AccumulatorHistory(new HistoryCollapsingDiscard(102, 3));
        frac2History.setTimeDataSource(stepCounter);
        DataPumpListenerSmart frac2Pump = new DataPumpListenerSmart(frac2DataSource, frac2History, dataInterval, (HistoryCollapsingDiscard)frac2History.getHistory());
        sim.integratorOS.getEventManager().addListener(frac2Pump);
        frac2History.setDataSink(fracPlot.getDataSet().makeDataSink());
        fracPlot.setLegend(new DataTag[]{frac2History.getTag()}, "optimal");
        
        DataSourceChiSlope slopeDataSource = new DataSourceChiSlope(dsvo);
        AccumulatorHistory slopeHistory = new AccumulatorHistory(new HistoryCollapsingDiscard(102, 3));
        slopeHistory.setTimeDataSource(stepCounter);
        DataPumpListenerSmart slopePump = new DataPumpListenerSmart(slopeDataSource, slopeHistory, dataInterval, (HistoryCollapsingDiscard)slopeHistory.getHistory());
        sim.integratorOS.getEventManager().addListener(slopePump);
        DisplayPlot slopePlot = new DisplayPlot();
        slopeHistory.setDataSink(slopePlot.getDataSet().makeDataSink());
        slopePlot.setLegend(new DataTag[]{slopeHistory.getTag()}, "slope");

        final DataSourceUa uaDataSource = new DataSourceUa(dsvo);
        AccumulatorHistory uaHistory = new AccumulatorHistory(new HistoryCollapsingDiscard(102, 3));
        uaHistory.setTimeDataSource(stepCounter);
        DataPumpListenerSmart uaPump = new DataPumpListenerSmart(uaDataSource, uaHistory, 10*dataInterval, (HistoryCollapsingDiscard)uaHistory.getHistory());
        sim.integratorOS.getEventManager().addListener(uaPump);
        DisplayPlot uaPlot = new DisplayPlot();
        uaHistory.setDataSink(uaPlot.getDataSet().makeDataSink());
        uaPlot.setLegend(new DataTag[]{uaHistory.getTag()}, "Ua");
        
        DataProcessorFunction dpRefUa2 = new DataProcessorFunction(new IFunction() {
            public double f(double x) {return x*x;}
        });
        forkA.addDataSink(dpRefUa2);
        final AccumulatorAverageFixed accRefUa2 = new AccumulatorAverageFixed();
        dpRefUa2.setDataSink(accRefUa2);

        DataProcessorFunction dpTargetUa2 = new DataProcessorFunction(new IFunction() {
            public double f(double x) {return x*x;}
        });
        forkB.addDataSink(dpTargetUa2);
        final AccumulatorAverageFixed accTargetUa2 = new AccumulatorAverageFixed();
        dpTargetUa2.setDataSink(accTargetUa2);
        DataSourceUa2 ua2DataSource = new DataSourceUa2(sim.integratorOS, dsvo, accRefUa2, accTargetUa2);
        AccumulatorHistory ua2History = new AccumulatorHistory(new HistoryCollapsingDiscard(102, 3));
        ua2History.setTimeDataSource(stepCounter);
        DataFork ua2Fork = new DataFork();
        DataPumpListenerSmart ua2Pump = new DataPumpListenerSmart(ua2DataSource, ua2Fork, 10*dataInterval, (HistoryCollapsingDiscard)ua2History.getHistory());
        ua2Fork.addDataSink(ua2History);
        sim.integratorOS.getEventManager().addListener(ua2Pump);
        ua2History.setPushInterval(1);
        ua2History.setDataSink(uaPlot.getDataSet().makeDataSink());
        uaPlot.setLegend(new DataTag[]{ua2History.getTag()}, "Ua2");

        final AccumulatorHistory aHistory = new AccumulatorHistory(new HistoryCollapsingDiscard(102, 3));
        DataProcessor aDataSource = new DataProcessAHT(uaDataSource, (HistoryCollapsingDiscard)aHistory.getHistory());
        ua2Fork.addDataSink(aDataSource);
        aHistory.setTimeDataSource(stepCounter);
        aDataSource.setDataSink(aHistory);
        DisplayPlot aPlot = new DisplayPlot();
        aHistory.setDataSink(slopePlot.getDataSet().makeDataSink());
        aPlot.setDoLegend(false);
        slopePlot.setLegend(new DataTag[]{aHistory.getTag()}, "a");
        
        
        DisplayPlot alphaChiPlot = new DisplayPlot();
        alphaChiPlot.getPlot().setXLog(true);
        alphaChiPlot.getPlot().setYLog(true);
        alphaChiPlot.setSize(350, 250);
        IDataSource alphaChi = new DataSourceAlphaChi(dsvo, 0);
        DataPumpListenerPow alphaChiPump = new DataPumpListenerPow(alphaChi, alphaChiPlot.getDataSet().makeDataSink(), 10*dataInterval, 10, 20);
        sim.integratorOS.getEventManager().addListener(alphaChiPump);
        IDataSource alphaAlpha = new DataSourceAlphaAlpha(dsvo);
        DataPumpListenerPow alphaAlphaPump = new DataPumpListenerPow(alphaAlpha, alphaChiPlot.getDataSet().makeDataSink(), 10*dataInterval, 10, 20);
        sim.integratorOS.getEventManager().addListener(alphaAlphaPump);

        DataSourceAlphaChi alphaChiErrP1 = new DataSourceAlphaChi(dsvo, 1);
        alphaChiErrP1.setSpiffiness(1);
        DataPumpListenerPow alphaChiPumpErrP1 = new DataPumpListenerPow(alphaChiErrP1, alphaChiPlot.getDataSet().makeDataSink(), 10*dataInterval, 10, 20);
        sim.integratorOS.getEventManager().addListener(alphaChiPumpErrP1);
        alphaChiPlot.setLegend(new DataTag[]{alphaChiErrP1.getTag()}, "chi+");
        
        DeviceSlider x0Slider = new DeviceSlider(sim.getController());
        final DeviceSlider omegaASlider = new DeviceSlider(sim.getController());
        final DeviceSlider omegaBSlider = new DeviceSlider(sim.getController());
        x0Slider.setShowValues(true);
        omegaASlider.setShowValues(true);
        omegaBSlider.setShowValues(true);
        x0Slider.setPrecision(1);
        omegaASlider.setPrecision(1);
        omegaBSlider.setPrecision(1);
        x0Slider.setEditValues(true);
        omegaASlider.setEditValues(true);
        omegaBSlider.setEditValues(true);
        Modifier x0Modifier = new Modifier() {
            public void setValue(double value) {
                sim.potentialB.setX0(new Vector1D(value));
            }
            public double getValue() {
                return sim.potentialB.getX0().getX(0);
            }
            public String getLabel() {return "x0";}
            public Dimension getDimension() {return Length.DIMENSION;}
        };
        x0Slider.setModifier(x0Modifier);
        x0Slider.setMinimum(0.0);
        x0Slider.setMaximum(3.0);
        x0Slider.setValue(2.0);
        omegaASlider.setModifier(new ModifierGeneral(sim.potentialA, "springConstant"));
        omegaASlider.setMinimum(0.1);
        omegaASlider.setMaximum(50.0);
        omegaASlider.setValue(10.0);
        omegaBSlider.setModifier(new ModifierGeneral(sim.potentialB, "springConstant"));
        omegaBSlider.setMinimum(0.1);
        omegaBSlider.setMaximum(10.0);
        omegaBSlider.setValue(1.0);
        
        final DeviceBox alphaCenterBox = new DeviceBox();
        alphaCenterBox.setController(sim.getController());
        alphaCenterBox.setModifier(new Modifier() {
            
            public void setValue(double newValue) {
                if (newValue <= 0) throw new IllegalArgumentException();
                double span = meterOverlapA.getAlphaSpan();
                meterOverlapA.setAlphaRange(newValue, span);
                meterOverlapB.setAlphaRange(newValue, span);
                for (int i=0; i<na; i++) {
                    ((AccumulatorAverageCollapsingLogAB)splitterB.getDataSink(i)).setMaxSample(1.0/meterOverlapB.getAlpha(i));
                }
            }
            
            public double getValue() {
                return meterOverlapA.getAlphaCenter();
            }
            
            public String getLabel() {
                return "alpha center";
            }
            
            public Dimension getDimension() {
                return Null.DIMENSION;
            }
        });

        final DataSourceScalar delta = new DataSourceScalar("exact",Energy.DIMENSION) {
            public double getDataAsScalar() {
                return 0.5*sim.boxA.getLeafList().size() * Math.log(omegaBSlider.getValue()/omegaASlider.getValue());
            }
        };
        
        
        DeviceBox alphaSpanBox = new DeviceBox();
        alphaSpanBox.setController(sim.getController());
        alphaSpanBox.setModifier(new Modifier() {
            
            public void setValue(double newValue) {
                if (newValue <= 0 || newValue > 100) throw new IllegalArgumentException();
                double center = meterOverlapA.getAlphaCenter();
                meterOverlapA.setAlphaRange(center, newValue);
                meterOverlapB.setAlphaRange(center, newValue);
                for (int i=0; i<na; i++) {
                    ((AccumulatorAverageCollapsingLogAB)splitterB.getDataSink(i)).setMaxSample(1.0/meterOverlapB.getAlpha(i));
                }

            }
            
            public double getValue() {
                return meterOverlapA.getAlphaSpan();
            }
            
            public String getLabel() {
                return "alpha span";
            }
            
            public Dimension getDimension() {
                return Null.DIMENSION;
            }
        });
        alphaSpanBox.setPrecision(6);

        DisplayPlot feRefTargetPlot = null;
        DisplayPlot stdevRefTargetPlot = new DisplayPlot();
        stdevRefTargetPlot.getPlot().setXLog(true);
        stdevRefTargetPlot.getPlot().setYLog(true);
        DisplayPlot biasRefTargetPlot = null;
        if (overlapAvgA instanceof DataSourceOverlapLogAvg) {
            dataSourceRefFE = new DataSourceAlphaFE(meterOverlapA, DataSourceAlphaFE.AVG);
            dataSourceRefFE.setOverlapSplitter(splitterA);
            dataSourceTargetFE = new DataSourceAlphaFE(meterOverlapA, DataSourceAlphaFE.AVG);
            dataSourceTargetFE.setOverlapSplitter(splitterB);
            feRefTargetPlot = new DisplayPlot();
            feRefTargetPlot.getPlot().setXLog(true);
            DataPumpListenerPow pump = new DataPumpListenerPow(dataSourceRefFE, feRefTargetPlot.getDataSet().makeDataSink(), 10*dataInterval, 10, 20);
            sim.integratorOS.getEventManager().addListener(pump);
            feRefTargetPlot.setLegend(new DataTag[]{dataSourceRefFE.getTag()}, "A");
            pump = new DataPumpListenerPow(dataSourceTargetFE, feRefTargetPlot.getDataSet().makeDataSink(), 10*dataInterval, 10, 20);
            sim.integratorOS.getEventManager().addListener(pump);
            feRefTargetPlot.setLegend(new DataTag[]{dataSourceTargetFE.getTag()}, "B");

            dataSourceRefStdev = new DataSourceAlphaFE(meterOverlapA, DataSourceAlphaFE.STDEV);
            dataSourceRefStdev.setOverlapSplitter(splitterA);
            dataSourceTargetStdev = new DataSourceAlphaFE(meterOverlapA, DataSourceAlphaFE.STDEV);
            dataSourceTargetStdev.setOverlapSplitter(splitterB);
            pump = new DataPumpListenerPow(dataSourceRefStdev, stdevRefTargetPlot.getDataSet().makeDataSink(), 10*dataInterval, 10, 20);
            sim.integratorOS.getEventManager().addListener(pump);
            stdevRefTargetPlot.setLegend(new DataTag[]{dataSourceRefStdev.getTag()}, "A");
            pump = new DataPumpListenerPow(dataSourceTargetStdev, stdevRefTargetPlot.getDataSet().makeDataSink(), 10*dataInterval, 10, 20);
            sim.integratorOS.getEventManager().addListener(pump);
            stdevRefTargetPlot.setLegend(new DataTag[]{dataSourceTargetStdev.getTag()}, "B");

            biasRefTargetPlot = new DisplayPlot();
            biasRefTargetPlot.getPlot().setXLog(true);
            biasRefTargetPlot.getPlot().setYLog(true);
            dataSourceRefBias = new DataSourceAlphaFE(meterOverlapA, DataSourceAlphaFE.BIAS);
            dataSourceRefBias.setOverlapSplitter(splitterA);
            dataSourceTargetBias = new DataSourceAlphaFE(meterOverlapA, DataSourceAlphaFE.BIAS);
            dataSourceTargetBias.setOverlapSplitter(splitterB);
            pump = new DataPumpListenerPow(dataSourceRefBias, biasRefTargetPlot.getDataSet().makeDataSink(), 10*dataInterval, 10, 20);
            sim.integratorOS.getEventManager().addListener(pump);
            biasRefTargetPlot.setLegend(new DataTag[]{dataSourceRefBias.getTag()}, "A");
            pump = new DataPumpListenerPow(dataSourceTargetBias, biasRefTargetPlot.getDataSet().makeDataSink(), 10*dataInterval, 10, 20);
            sim.integratorOS.getEventManager().addListener(pump);
            biasRefTargetPlot.setLegend(new DataTag[]{dataSourceTargetBias.getTag()}, "B");
        }
        
        DataSourceScalar refAntibiasDataSource = new DataSourceScalar("antibias", Null.DIMENSION) {
            public double getDataAsScalar() {
                double refChi = dsvo.getAverageAndError(true, alphaChoice)[0];
                long refCount = sim.integratorA.getStepCount();
                double antibiasRefChi = (refChi * refCount + 1) / (refCount+1);
                return Math.log(antibiasRefChi);
            }
        };
        AccumulatorHistory refAntibiasHistory = new AccumulatorHistory(new HistoryComplete());
        refAntibiasHistory.setTimeDataSource(new DataSourceCountSteps(sim.integratorA));
        DataPumpListenerPowStrict refAntibiasPump = new DataPumpListenerPowStrict(refAntibiasDataSource, refAntibiasHistory);
        sim.integratorA.getEventManager().addListener(refAntibiasPump);
        refAntibiasHistory.setPushInterval(1);
        refAntibiasHistory.addDataSink(feRefTargetPlot.getDataSet().makeDataSink());
        feRefTargetPlot.setLegend(new DataTag[]{refAntibiasHistory.getTag()}, "A(ab)");

        DataSourceScalar targetAntibiasDataSource = new DataSourceScalar("antibias", Null.DIMENSION) {
            public double getDataAsScalar() {
                double targetChi = dsvo.getAverageAndError(false, alphaChoice)[0];
                long targetCount = sim.integratorB.getStepCount();
                double antibiasRefChi = (targetChi * targetCount + 1.0/alphaChoice) / (targetCount+1);
                return Math.log(antibiasRefChi);
            }
        };
        AccumulatorHistory targetAntibiasHistory = new AccumulatorHistory(new HistoryComplete());
        targetAntibiasHistory.setTimeDataSource(new DataSourceCountSteps(sim.integratorB));
        DataPumpListenerPowStrict targetAntibiasPump = new DataPumpListenerPowStrict(targetAntibiasDataSource, targetAntibiasHistory);
        sim.integratorB.getEventManager().addListener(targetAntibiasPump);
        targetAntibiasHistory.setPushInterval(1);
        targetAntibiasHistory.addDataSink(feRefTargetPlot.getDataSet().makeDataSink());
        feRefTargetPlot.setLegend(new DataTag[]{targetAntibiasHistory.getTag()}, "B(ab)");
        

        sim.integratorA.getEventManager().addListener(new IntegratorListener() {
            public void integratorStepStarted(IntegratorEvent e) {}
            public void integratorInitialized(IntegratorEvent e) {}
            
            public void integratorStepFinished(IntegratorEvent e) {
                long stepCount = sim.integratorA.getStepCount();
                int ls = (int)Math.round(Math.log(stepCount)/Math.log(2));
                if (ls < 10 || stepCount != 1L<<ls) return;
                int numAlphas = meterOverlapA.getNumAlpha();
                double[][] avgs = new double[numAlphas][ls+1];
                double[][] sdevs = new double[numAlphas][ls+1];
                double[][] abavgs = new double[numAlphas][ls+1];
                double[][] absdevs = new double[numAlphas][ls+1];
                for (int i=0; i<numAlphas; i++) {
                    IData davg = ((AccumulatorAverageCollapsingLog)splitterA.getDataSink(i)).getAverageLogs();
                    IData dsdev = ((AccumulatorAverageCollapsingLog)splitterA.getDataSink(i)).getStdevLog();
                    for (int j=0; j<=ls; j++) {
                        avgs[i][j] = davg.getValue(j);
                        sdevs[i][j] = dsdev.getValue(j);
                    }
                    IData dabavg = ((AccumulatorAverageCollapsingLogAB)splitterA.getDataSink(i)).getAverageAntibiasedLogs();
                    IData dabsdev = ((AccumulatorAverageCollapsingLogAB)splitterA.getDataSink(i)).getStdevantibiasedLog();
                    for (int j=0; j<=ls; j++) {
                        abavgs[i][j] = dabavg.getValue(j);
                        absdevs[i][j] = dabsdev.getValue(j);
                    }
                }
                for (int j=0; j<=ls; j++) {
                    long steps = 1L<<j;
                    System.out.print("A "+steps);
                    for (int i=0; i<numAlphas; i++) {
                        System.out.print("   "+avgs[i][j]+" "+sdevs[i][j]+" "+abavgs[i][j]+" "+absdevs[i][j]);
                    }
                    System.out.println();
                }
            }
        });
        sim.integratorB.getEventManager().addListener(new IntegratorListener() {
            public void integratorStepStarted(IntegratorEvent e) {}
            public void integratorInitialized(IntegratorEvent e) {}
            
            public void integratorStepFinished(IntegratorEvent e) {
                long stepCount = sim.integratorB.getStepCount();
                int ls = (int)Math.round(Math.log(stepCount)/Math.log(2));
                if (ls < 10 || stepCount != 1L<<ls) return;
                int numAlphas = meterOverlapA.getNumAlpha();
                double[][] avgs = new double[numAlphas][ls+1];
                double[][] sdevs = new double[numAlphas][ls+1];
                double[][] abavgs = new double[numAlphas][ls+1];
                double[][] absdevs = new double[numAlphas][ls+1];
                for (int i=0; i<numAlphas; i++) {
                    IData davg = ((AccumulatorAverageCollapsingLog)splitterB.getDataSink(i)).getAverageLogs();
                    IData dsdev = ((AccumulatorAverageCollapsingLog)splitterB.getDataSink(i)).getStdevLog();
                    for (int j=0; j<=ls; j++) {
                        avgs[i][j] = davg.getValue(j);
                        sdevs[i][j] = dsdev.getValue(j);
                    }
                    IData dabavg = ((AccumulatorAverageCollapsingLogAB)splitterB.getDataSink(i)).getAverageAntibiasedLogs();
                    IData dabsdev = ((AccumulatorAverageCollapsingLogAB)splitterB.getDataSink(i)).getStdevantibiasedLog();
                    for (int j=0; j<=ls; j++) {
                        abavgs[i][j] = dabavg.getValue(j);
                        absdevs[i][j] = dabsdev.getValue(j);
                    }
                }
                for (int j=0; j<=ls; j++) {
                    long steps = 1L<<j;
                    System.out.print("B "+steps);
                    for (int i=0; i<numAlphas; i++) {
                        System.out.print("   "+avgs[i][j]+" "+sdevs[i][j]+" "+abavgs[i][j]+" "+absdevs[i][j]);
                    }
                    System.out.println();
                }
            }
        });

        DataSourceScalar refErrDataSource = new DataSourceScalar("A", Null.DIMENSION) {
            public double getDataAsScalar() {
                double err = dsvo.getLogAverageAndError(true, alphaChoice)[1];
                return err > 0 ? err : Double.NaN;
            }
        };
        AccumulatorHistory refErrHistory = new AccumulatorHistory(new HistoryComplete());
        refErrHistory.setTimeDataSource(new DataSourceCountSteps(sim.integratorA));
        DataPumpListenerPowStrict refErrPump = new DataPumpListenerPowStrict(refErrDataSource, refErrHistory);
        sim.integratorA.getEventManager().addListener(refErrPump);
        refErrHistory.setPushInterval(1);
        refErrHistory.addDataSink(stdevRefTargetPlot.getDataSet().makeDataSink());
        stdevRefTargetPlot.setLegend(new DataTag[]{refErrHistory.getTag()}, "A");

        DataSourceScalar targetErrDataSource = new DataSourceScalar("B", Null.DIMENSION) {
            public double getDataAsScalar() {
                double err = dsvo.getLogAverageAndError(false, alphaChoice)[1];
                return err > 0 ? err : Double.NaN;
            }
        };
        AccumulatorHistory targetErrHistory = new AccumulatorHistory(new HistoryComplete());
        DataPumpListenerPowStrict targetErrPump = new DataPumpListenerPowStrict(targetErrDataSource, targetErrHistory);
        sim.integratorB.getEventManager().addListener(targetErrPump);
        targetErrHistory.setTimeDataSource(new DataSourceCountSteps(sim.integratorB));
        targetErrHistory.setPushInterval(1);
        targetErrHistory.addDataSink(stdevRefTargetPlot.getDataSet().makeDataSink());
        stdevRefTargetPlot.setLegend(new DataTag[]{targetErrHistory.getTag()}, "B");

//        sim.integratorOS.getEventManager().addListener(new IntegratorListener() {
//            public void integratorStepStarted(IIntegratorEvent e) {}
//            public void integratorStepFinished(IIntegratorEvent e) {
//                if (sim.integratorA.getStepCount() > 1<<22 && sim.integratorB.getStepCount() > 1<<22) System.exit(0);
//            }
//            public void integratorInitialized(IIntegratorEvent e) {}
//        });
        
        DeviceBox numAlphaBox = new DeviceBox();
        numAlphaBox.setController(sim.getController());
        numAlphaBox.setInteger(true);
        numAlphaBox.setModifier(new Modifier() {
            
            public void setValue(double newValue) {
                if (newValue <= 0 || newValue >= 50) throw new IllegalArgumentException();
                int numAlpha = (int)newValue;
                meterOverlapA.setNumAlpha(numAlpha);
                meterOverlapB.setNumAlpha(numAlpha);
                forkA.putDataInfo(meterOverlapA.getDataInfo());
                for (int i=0; i<numAlpha; i++) {
                    ((AccumulatorAverageCollapsingLogAB)splitterA.getDataSink(i)).setMaxSample(1);
                }
                forkB.putDataInfo(meterOverlapB.getDataInfo());
                for (int i=0; i<numAlpha; i++) {
                    ((AccumulatorAverageCollapsingLogAB)splitterB.getDataSink(i)).setMaxSample(1.0/meterOverlapB.getAlpha(i));
                }
                for (DataPumpListener pump : allPumps) {
                    sim.integratorOS.getEventManager().removeListener(pump);
                }
                allPumps.clear();
                fePlot.getDataSet().reset();
            }
            
            public double getValue() {
                return meterOverlapA.getNumAlpha();
            }
            
            public String getLabel() {
                return "# alpha";
            }
            
            public Dimension getDimension() {
                return Null.DIMENSION;
            }
        });

        final DeviceBox alphaChoiceBox = new DeviceBox();
        alphaChoiceBox.setModifier(new Modifier() {
            public void setValue(double newValue) {
                AlphaSource alphaSource = dsvo.getAlphaSource();
                int numAlpha = alphaSource.getNumAlpha();
                if (newValue < alphaSource.getAlpha(0) || newValue > alphaSource.getAlpha(numAlpha-1)) {
                    throw new IllegalArgumentException();
                }
                alphaChoice = newValue;
                dataSourceRefFE.setAlpha(alphaChoice);
                dataSourceTargetFE.setAlpha(alphaChoice);
                dataSourceRefStdev.setAlpha(alphaChoice);
                dataSourceTargetStdev.setAlpha(alphaChoice);
                dataSourceRefBias.setAlpha(alphaChoice);
                dataSourceTargetBias.setAlpha(alphaChoice);
            }
            
            public double getValue() {
                return alphaChoice;
            }
            
            public String getLabel() {
                return "alpha";
            }
            
            public Dimension getDimension() {
                return Null.DIMENSION;
            }
        });


        final DeviceNSelector nSlider = new DeviceNSelector(sim.getController());
        nSlider.setBox(sim.boxA);
        nSlider.setSpecies(sim.species);
        nSlider.setMaximum(50);
        nSlider.setMinimum(1);
        nSlider.setShowValues(true);
        nSlider.setLabel("Number of atoms");
        nSlider.setShowBorder(true);

        AccumulatorHistory deltaHistory = new AccumulatorHistory(new HistoryCollapsingDiscard(102, 3));
        DataPumpListenerSmart exactPump = new DataPumpListenerSmart(delta, deltaHistory, dataInterval, (HistoryCollapsingDiscard)deltaHistory.getHistory());
        deltaHistory.setDataSink(fePlot.getDataSet().makeDataSink());
        sim.integratorOS.getEventManager().addListener(exactPump);
        dataStreamPumps.add(exactPump);
        deltaHistory.setTimeDataSource(stepCounter);

        AccumulatorHistory deltaHistoryLn = new AccumulatorHistory(new HistoryComplete());
        DataPumpListenerPowStrict exactPumpLn = new DataPumpListenerPowStrict(delta, deltaHistoryLn);
        deltaHistoryLn.setDataSink(feLnPlot.getDataSet().makeDataSink());
        sim.integratorOS.getEventManager().addListener(exactPumpLn);
        dataStreamPumps.add(exactPumpLn);
        deltaHistoryLn.setTimeDataSource(stepCounter);

        final DisplayPlot uPlot = new DisplayPlot();
        final double yMax = 2.0;
        uPlot.getPlot().setYRange(0.0, yMax);
        
        Function fUA = new Function() {
            public double f(double x) {
                double x0 = sim.potentialA.getX0().getX(0);
                return 0.5*sim.potentialA.getSpringConstant()*(x - x0)*(x - x0);
            }
        };
        Function fUB = new Function() {
            public double f(double x) {
                double x0 = sim.potentialB.getX0().getX(0);
                return 0.5*sim.potentialB.getSpringConstant()*(x - x0)*(x - x0);
            }
        };

        final DataSourceFunction uA = new DataSourceFunction("A",Null.DIMENSION,fUA,100,"x",Length.DIMENSION);
        final DataSourceFunction uB = new DataSourceFunction("B",Null.DIMENSION,fUB,100,"x",Length.DIMENSION);
        double Lx = sim.boxA.getBoundary().getBoxSize().getX(0);
        uA.getXSource().setXMin(-Lx);
        uB.getXSource().setXMin(-Lx);
        uA.getXSource().setXMax(Lx);
        uB.getXSource().setXMax(Lx);
        final DataPump uAPump = new DataPump(uA, uPlot.getDataSet().makeDataSink());
        final DataPump uBPump = new DataPump(uB, uPlot.getDataSet().makeDataSink());
        final IAction uUpdate = new IAction() {
            public void actionPerformed() {
                uA.update();
                uB.update();
                uAPump.actionPerformed();
                uBPump.actionPerformed();
                
                double span = meterOverlapA.getAlphaSpan();
                double alpha = Math.exp(-delta.getDataAsScalar());
                meterOverlapA.setAlphaRange(alpha, span);
                meterOverlapB.setAlphaRange(alpha, span);
                for (int i=0; i<na; i++) {
                    ((AccumulatorAverageCollapsingLogAB)splitterB.getDataSink(i)).setMaxSample(1.0/meterOverlapB.getAlpha(i));
                }
                alphaCenterBox.doUpdate();
                alphaChoice = alpha;
                alphaChoiceBox.doUpdate();

                dataSourceRefFE.setAlpha(alphaChoice);
                dataSourceTargetFE.setAlpha(alphaChoice);
                dataSourceRefStdev.setAlpha(alphaChoice);
                dataSourceTargetStdev.setAlpha(alphaChoice);
                dataSourceRefBias.setAlpha(alphaChoice);
                dataSourceTargetBias.setAlpha(alphaChoice);
            }
        };
        omegaASlider.setPostAction(uUpdate);
        omegaBSlider.setPostAction(uUpdate);
        x0Slider.setPostAction(uUpdate);
        nSlider.setPostAction(new IAction() {
            public void actionPerformed() {
                sim.boxB.setNMolecules(sim.species, (int)nSlider.getValue());
                uUpdate.actionPerformed();
            }
        });

        uPlot.getDataSet().setUpdatingOnAnyChange(true);

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();
        GridBagConstraints horizGBC = SimulationPanel.getHorizGBC();

        //controls -- start/pause and sliders
        JTabbedPane sliderPanel = new JTabbedPane();
        sliderPanel.add(x0Slider.graphic(), "x0");
        sliderPanel.add(omegaASlider.graphic(), "omegaA");
        sliderPanel.add(omegaBSlider.graphic(), "omegaB");
        
        JPanel topPanel = new JPanel(new GridBagLayout());
        JPanel displayPanel = new JPanel(new GridBagLayout());
        topPanel.add(sliderPanel, vertGBC);
        topPanel.add(nSlider.graphic(), vertGBC);
        JPanel alphaPanel = new JPanel(new GridLayout(2, 2));
        alphaPanel.add(alphaCenterBox.graphic());
        alphaPanel.add(alphaSpanBox.graphic());
        alphaPanel.add(numAlphaBox.graphic());
        alphaPanel.add(alphaChoiceBox.graphic());
        topPanel.add(alphaPanel, vertGBC);
        displayPanel.add(uPlot.graphic(), vertGBC);
        JTabbedPane displayBoxPanel = new JTabbedPane();
        displayBoxPanel.add(displayBoxA.graphic(), "box A");
        displayBoxPanel.add(displayBoxB.graphic(), "box B");
        displayPanel.add(displayBoxPanel, vertGBC);
        GridBagConstraints mygbc = new GridBagConstraints();
        mygbc.gridx = 1;
        mygbc.gridy = 0;
        mygbc.gridheight = 3;
        topPanel.add(displayPanel, mygbc);

        stdevRefTargetPlot.setSize(350,250);
        if (feRefTargetPlot != null) {
            feRefTargetPlot.setSize(350,250);
            biasRefTargetPlot.setSize(350,250);
        }
        aPlot.setSize(450,250);
        uaPlot.setSize(450,250);
        slopePlot.setSize(450,250);
        fracPlot.setSize(450,250);
        fePlot.setSize(450, 250);
        feLnPlot.setSize(450, 250);
        uPlot.setSize(450, 250);
        if (dWPlot != null) {
            dWPlot.setSize(450, 250);
        }

        getPanel().graphicsPanel.add(topPanel);

        JTabbedPane plotTabs = new JTabbedPane();

        getPanel().graphicsPanel.add(plotTabs, vertGBC);

        getController().getReinitButton().setPostAction(new IAction() {
        	public void actionPerformed() {
                displayBoxA.repaint();
                displayBoxB.repaint();
                fePlot.getPlot().repaint();
                feLnPlot.getPlot().repaint();
        	}
        });

        uUpdate.actionPerformed();

        JPanel tab2 = new JPanel(new GridBagLayout());
        plotTabs.add(tab2, "Free Energy");
        tab2.add(fePlot.graphic(), horizGBC);
        JPanel tab2a = new JPanel(new GridBagLayout());
        plotTabs.add(tab2a, "Free Energy (log)");
        tab2a.add(feLnPlot.graphic(), horizGBC);
        JPanel tab3 = new JPanel(new GridBagLayout());
        plotTabs.add(tab3, "Alpha vs. Chi");
        tab3.add(alphaChiPlot.graphic(), horizGBC);
        JPanel tab4 = new JPanel(new GridBagLayout());
        plotTabs.add(tab4, "Frac A");
        tab4.add(fracPlot.graphic(), horizGBC);
        JPanel tab5 = new JPanel(new GridBagLayout());
        plotTabs.add(tab5, "Chi Slope");
        tab5.add(slopePlot.graphic(), horizGBC);
        JPanel tab6 = new JPanel(new GridBagLayout());
        plotTabs.add(tab6, "Ua");
        tab6.add(uaPlot.graphic(), horizGBC);
        if (feRefTargetPlot != null) {
            JPanel tab7 = new JPanel(new GridBagLayout());
            plotTabs.add(tab7, "AB FE");
            tab7.add(feRefTargetPlot.graphic(), horizGBC);
        }
        JPanel tab8 = new JPanel(new GridBagLayout());
        plotTabs.add(tab8, "AB Stdev");
        tab8.add(stdevRefTargetPlot.graphic(), horizGBC);
        if (biasRefTargetPlot != null) {
            JPanel tab10 = new JPanel(new GridBagLayout());
            plotTabs.add(tab10, "AB bias");
            tab10.add(biasRefTargetPlot.graphic(), horizGBC);
        }
        
        if (dWPlot != null) {
            JPanel tab13 = new JPanel(new GridBagLayout());
            plotTabs.add(tab13, "dW");
            tab13.add(dWPlot.graphic(), horizGBC);
        }
        
        final DataSinkReweightedDeltaU deltaUA = new DataSinkReweightedDeltaU(meterOverlapA);
        forkA.addDataSink(deltaUA);
        DisplayTextBox sWABox = new DisplayTextBox();
        DataPumpListenerPow sWAPump = new DataPumpListenerPow(new DataSourceScalar("sWA", Null.DIMENSION) {
            public double getDataAsScalar() {
                // sWA = <UA - UW>W - (FA-FW)
                return deltaUA.getDeltaU(alphaChoice) - Math.log(dsvo.getAverageAndError(true, alphaChoice)[0]);
            }
        }, sWABox, dataInterval*1000, 10, 20);
        sim.integratorA.getEventManager().addListener(sWAPump);

        final DataSinkReweightedDeltaU deltaUB = new DataSinkReweightedDeltaU(meterOverlapB);
        forkB.addDataSink(deltaUB);
        DisplayTextBox sWBBox = new DisplayTextBox();
        DataPumpListenerPow sWBPump = new DataPumpListenerPow(new DataSourceScalar("sWB", Null.DIMENSION) {
            public double getDataAsScalar() {
                // sWB = <UB - UW>W - (FB-FW)
                return deltaUB.getDeltaU(alphaChoice) - Math.log(dsvo.getAverageAndError(false, alphaChoice)[0]);
            }
        }, sWBBox, dataInterval*1000, 10, 20);
        sim.integratorB.getEventManager().addListener(sWBPump);

        DisplayTextBox MABox = new DisplayTextBox();
        DataSourceScalar dataSourceMA = new DataSourceScalar("MA", Quantity.DIMENSION) {
            public double getDataAsScalar() {
                double sWA = deltaUA.getDeltaU(alphaChoice) - Math.log(dsvo.getAverageAndError(true, alphaChoice)[0]);
                return 1.0 + Math.exp(sWA)*Math.sqrt(2*Math.PI);
            }
        };
        DataPumpListenerPow pumpMA = new DataPumpListenerPow(dataSourceMA, MABox, dataInterval*1000, 10, 20);
        sim.integratorA.getEventManager().addListener(pumpMA);
        AccumulatorHistory historyMA = new AccumulatorHistory(new HistoryComplete());
        DataPumpListenerPowStrict pumpMALog = new DataPumpListenerPowStrict(dataSourceMA, historyMA);
        historyMA.setTimeDataSource(new DataSourceCountSteps(sim.integratorA));
        sim.integratorA.getEventManager().addListener(pumpMALog);

        DisplayTextBox MBBox = new DisplayTextBox();
        DataSourceScalar dataSourceMB = new DataSourceScalar("MB", Quantity.DIMENSION) {
            public double getDataAsScalar() {
                double sWB = deltaUB.getDeltaU(alphaChoice) - Math.log(dsvo.getAverageAndError(false, alphaChoice)[0]);
                return 1.0 + Math.exp(sWB)*Math.sqrt(2*Math.PI);
            }
        };
        DataPumpListenerPow pumpMB = new DataPumpListenerPow(dataSourceMB, MBBox, dataInterval*1000, 10, 20);
        sim.integratorB.getEventManager().addListener(pumpMB);
        AccumulatorHistory historyMB = new AccumulatorHistory(new HistoryComplete());
        DataPumpListenerPowStrict pumpMBLog = new DataPumpListenerPowStrict(dataSourceMB, historyMB);
        historyMB.setTimeDataSource(new DataSourceCountSteps(sim.integratorB));
        sim.integratorB.getEventManager().addListener(pumpMBLog);

        
        JPanel tab11 = new JPanel(new GridBagLayout());
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        tab11.add(sWABox.graphic(), gbc);
        gbc.gridx = 1;
        tab11.add(sWBBox.graphic(), gbc);
        gbc.gridx = 0;
        gbc.gridy = 1;
        tab11.add(MABox.graphic(), gbc);
        gbc.gridx = 1;
        tab11.add(MBBox.graphic(), gbc);
        plotTabs.add(tab11, "s");
        

        DataSourceScalar dataSourceMM = new DataSourceScalar("MM", Quantity.DIMENSION) {
            public double getDataAsScalar() {
                return sim.integratorB.getStepCount() + sim.integratorA.getStepCount();
            }
        };
        AccumulatorHistory historyMM = new AccumulatorHistory(new HistoryComplete());
        DataPumpListenerPowStrict pumpMMLog = new DataPumpListenerPowStrict(dataSourceMM, historyMM);
        historyMM.setTimeDataSource(dataSourceMM);
        sim.integratorOS.getEventManager().addListener(pumpMMLog);

        DisplayPlot plotM = new DisplayPlot();
        historyMA.setDataSink(plotM.getDataSet().makeDataSink());
        historyMB.setDataSink(plotM.getDataSet().makeDataSink());
        historyMM.setDataSink(plotM.getDataSet().makeDataSink());
        plotM.setLegend(new DataTag[]{historyMA.getTag()}, "A");
        plotM.setLegend(new DataTag[]{historyMB.getTag()}, "B");
        plotM.setLegend(new DataTag[]{historyMM.getTag()}, "M");
        plotM.getPlot().setXLog(true);
        plotM.getPlot().setYLog(true);
        plotM.setSize(450,250);
        JPanel tab12 = new JPanel();
        tab12.add(plotM.graphic(), horizGBC);
        plotTabs.add(tab12, "M");
        
        sim.activityIntegrate.setSleepPeriod(0);
    }

    public static void main(String[] args) {
        final MultiharmonicMC sim = new MultiharmonicMC();
        MultiharmonicGraphicMC simGraphic = new MultiharmonicGraphicMC(sim, sim.getSpace());
        SimulationGraphic.makeAndDisplayFrame(simGraphic.getPanel(), APP_NAME);
    }

    public static class DataProcessAHT extends DataProcessor {
        private final DataSourceUa uaDataSource;
        private final HistoryCollapsingDiscard aHistory;
        protected final DataDouble data = new DataDouble();

        public DataProcessAHT(DataSourceUa uaDataSource,
                HistoryCollapsingDiscard aHistory) {
            this.uaDataSource = uaDataSource;
            this.aHistory = aHistory;
            dataInfo = new DataInfoDouble("a", Null.DIMENSION);
        }

        protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
            return dataInfo;
        }

        protected IData processData(IData inputData) {
            if (aHistory.willDiscardNextData()) return data;
            double ua2 = inputData.getValue(0);
            double ua = uaDataSource.getDataAsScalar();
            data.x = -(ua - ua2) / ua;
            return data;
        }
    }

    public static class DataSourceFunctionRaw implements IDataSource, DataSourceIndependent {
        protected final DataFunction data;
        protected final DataInfoFunction dataInfo;
        protected final DataDoubleArray xData;
        protected final DataInfoDoubleArray xDataInfo;
        protected final DataTag tag, xTag;
        
        public DataSourceFunctionRaw(double[] x, double[] y) {
            xData = new DataDoubleArray(new int[]{x.length}, x);
            xDataInfo = new DataInfoDoubleArray("x", Quantity.DIMENSION, new int[]{x.length});
            xTag = new DataTag();
            xDataInfo.addTag(xTag);
            data = new DataFunction(new int[]{y.length}, y);
            dataInfo = new DataInfoFunction("y", Null.DIMENSION, this);
            tag = new DataTag();
            dataInfo.addTag(tag);
        }

        public IData getData() {return data;}
        public DataDoubleArray getIndependentData(int i) {return xData;}
        public DataInfoDoubleArray getIndependentDataInfo(int i) {return xDataInfo;}
        public int getIndependentArrayDimension() {return 1;}
        public DataTag getIndependentTag() {return xTag;}
        public DataTag getTag() {return tag;}
        public IDataInfo getDataInfo() {return dataInfo;}
    }

    public static class DataSourceAlphaFE implements IDataSource, DataSourceIndependent {

        protected DataSplitter splitter;
        protected double[][] allData;
        protected double[] lnN;
        protected DataDoubleArray nData;
        protected DataInfoDoubleArray nDataInfo;
        protected DataFunction data;
        protected DataInfoFunction dataInfo;
        protected final AkimaSpline akima = new AkimaSpline();
        protected double[] alpha = new double[1];
        protected double[] allLnAlpha;
        protected final DataTag tag, nTag;
        protected final AlphaSource alphaSource;
        protected final int which;
        public static final int AVG = 1, STDEV = 2, SKEW = 4, BIAS = 5;

        public DataSourceAlphaFE(AlphaSource alphaSource, int which) {
            this.alphaSource = alphaSource;
            this.which = which;
            tag = new DataTag();
            nTag = new DataTag();
            nData = new DataDoubleArray(0);
            nDataInfo = new DataInfoDoubleArray("N", Quantity.DIMENSION, new int[]{0});
            nDataInfo.addTag(nTag);
            data = new DataFunction(new int[]{0});
            dataInfo = new DataInfoFunction("fe", Null.DIMENSION, this);
            dataInfo.addTag(tag);
            allData = new double[0][0];
            allLnAlpha = new double[0];
        }
        
        public void setOverlapSplitter(DataSplitter splitter) {
            this.splitter = splitter;
        }
        
        public void setAlpha(double newAlpha) {
            alpha[0] = Math.log(newAlpha); 
        }
        
        public IData getData() {
            AccumulatorAverageCollapsingLog[] acc = new AccumulatorAverageCollapsingLog[splitter.getNumDataSinks()];
            for (int i=0; i<acc.length; i++) {
                acc[i] = (AccumulatorAverageCollapsingLog)splitter.getDataSink(i);
            }
            if (allLnAlpha.length != acc.length) {
                allLnAlpha = new double[splitter.getNumDataSinks()];
                allData = new double[0][0];
            }
            int nDoubles = acc[0].getDataInfo().getLength();
            if (allData.length != nDoubles) {
                nData = new DataDoubleArray(nDoubles);
                nDataInfo = new DataInfoDoubleArray("N", Quantity.DIMENSION, new int[]{nDoubles});
                nDataInfo.addTag(nTag);
                data = new DataFunction(new int[]{nDoubles});
                dataInfo = new DataInfoFunction("fe", Null.DIMENSION, this);
                dataInfo.addTag(tag);

                allData = new double[nDoubles][0];
                for (int i=0; i<nDoubles; i++) {
                    allData[i] = new double[acc.length];
                }
                lnN = new double[nDoubles];
                double[] n = nData.getData();
                double ln2 = Math.log(2);
                for (int j=0; j<nDoubles; j++) {
                    lnN[j] = j*ln2;
                    n[j] = 1L<<j;
                }
            }
            for (int i=0; i<acc.length; i++) {
                IData accData = null;
                switch (which) {
                    case AVG:
                        accData = acc[i].getAverageLogs();
                        break;
                    case STDEV:
                        accData = acc[i].getStdevLog();
                        break;
                    case SKEW:
                        accData = acc[i].getSkewLog();
                        break;
                    case BIAS:
                        accData = acc[i].getAverageLogs();
                        for (int j=0; j<allData.length; j++) {
                            allData[j][i] = -accData.getValue(j);
                        }
                        accData = acc[i].getAverages();
                        for (int j=0; j<allData.length; j++) {
                            allData[j][i] += Math.log(accData.getValue(j));
                            if (allData[j][i] <= 0) allData[j][i] = Double.NaN; 
                        }
                        break;
                    default:
                        throw new RuntimeException("oops");
                }
                if (which != BIAS) {
                    for (int j=0; j<allData.length; j++) {
                        allData[j][i] = accData.getValue(j);
                    }
                }
            }

            double[] y = data.getData();
            for (int i=0; i<acc.length; i++) {
                allLnAlpha[i] = Math.log(alphaSource.getAlpha(i));
            }
            for (int i=0; i<y.length; i++) {
                if (acc.length == 1) {
                    y[i] = allData[i][0];
                }
                else {
                    akima.setInputData(allLnAlpha, allData[i]);
                    y[i] = akima.doInterpolation(alpha)[0];
                }
            }
            return data;
        }

        public DataTag getTag() {
            return tag;
        }

        public IDataInfo getDataInfo() {
            return dataInfo;
        }

        public DataDoubleArray getIndependentData(int i) {
            return nData;
        }

        public DataInfoDoubleArray getIndependentDataInfo(int i) {
            return nDataInfo;
        }

        public int getIndependentArrayDimension() {
            return 1;
        }

        public DataTag getIndependentTag() {
            return nTag;
        }
    }
    
    public static class DataSourceBias extends DataSourceScalar {

        protected AccumulatorAverageCollapsingLog[] acc;
        protected double[] allData;
        protected final AkimaSpline akima = new AkimaSpline();
        protected final AlphaSource alphaSource;
        protected final double[] alpha = new double[1];
        protected double[] allLnAlpha;
        public static final int AVG = 1, STDEV = 2, STDEV2 = 3, SKEW = 4, BIAS = 5;

        public DataSourceBias(AlphaSource alphaSource) {
            super("bias", Null.DIMENSION);
            this.alphaSource = alphaSource;
            allData = new double[0];
        }
        
        public void setAccumulators(AccumulatorAverageCollapsingLog[] newAcc) {
            acc = newAcc;
            allLnAlpha = new double[acc.length];
            allData = new double[acc.length];
        }
        
        public void setAlpha(double newAlpha) {
            alpha[0] = Math.log(newAlpha); 
        }
        
        public double getDataAsScalar() {
            int nDoubles = acc[0].getDataInfo().getLength();
            if (nDoubles < 2) return Double.NaN;
            for (int i=0; i<acc.length; i++) {
                allData[i] = -acc[i].getAverageLogs().getValue(nDoubles-1);
                allData[i] += Math.log(acc[i].getAverages().getValue(nDoubles-1));
            }

            for (int i=0; i<acc.length; i++) {
                allLnAlpha[i] = Math.log(alphaSource.getAlpha(i));
            }
            akima.setInputData(allLnAlpha, allData);
            double bias = akima.doInterpolation(alpha)[0];
            if (bias <= 0) return Double.NaN;
            return bias;
        }
    }
    
    public static class DataSourceChiSlope extends DataSourceScalar {
        protected final DataOverlap dsvo;
        
        public DataSourceChiSlope(DataOverlap dsvo) {
            super("slope", Null.DIMENSION);
            this.dsvo = dsvo;
        }

        public double getDataAsScalar() {
            int nBennetPoints = dsvo.getAlphaSource().getNumAlpha();
            if (nBennetPoints == 1) {
                return Double.NaN;
            }

            double[] lnAlpha = new double[nBennetPoints];
            double[] lnAlphaDiff = new double[nBennetPoints];

            for (int j=0; j<nBennetPoints; j++) {
                double refOverlap = dsvo.getRefSource().getAverage(j);
                double targetOverlap = dsvo.getTargetSource().getAverage(j);
                lnAlphaDiff[j] += Math.log(refOverlap/targetOverlap);

                double jAlpha = dsvo.getAlphaSource().getAlpha(j);
                lnAlpha[j] = Math.log(jAlpha);
                lnAlphaDiff[j] -= lnAlpha[j];
            }

            double slope = 0;
            if (lnAlphaDiff[0] < 0) {
                // first new alpha is less than initial first alpha
                return Double.NaN;
            }
            else if (lnAlphaDiff[nBennetPoints-1] > 0) {
                return Double.NaN;
            }
            else if (nBennetPoints > 4) {
                AkimaSpline spline = new AkimaSpline();
                spline.setInputData(lnAlpha, lnAlphaDiff);
                double min = lnAlpha[0];
                double max = lnAlpha[nBennetPoints-1];
                double ymin = lnAlphaDiff[0];
                double ymax = lnAlphaDiff[nBennetPoints-1];
                double[] x = new double[1];
                x[0] = min+(max-min)/(ymax-ymin)*(-ymin);
                while (x[0] > min && x[0] < max) {
                    x[0] = 0.5 * (min + max);
                    double y = spline.doInterpolation(x)[0];
                    if (y == 0) {
                        break;
                    }
                    if (y*min > 0) {
                        min = x[0];
                        ymin = y;
                    }
                    else {
                        max = x[0];
                        ymax = y;
                    }
                    x[0] = min+(max-min)/(ymax-ymin)*(-ymin);
                }
                slope = spline.doInterpolationDy(x)[0] + 1;
            }
            else {
                //linear interpolation (only 3 points)
                for (int i=0; i<nBennetPoints; i++) {
                    if (lnAlphaDiff[i] > 0 && lnAlphaDiff[i+1] < 0) {
                        slope = ((lnAlphaDiff[i+1] + lnAlpha[i+1]) - (lnAlphaDiff[i] + lnAlpha[i])) / (lnAlpha[i+1] - lnAlpha[i]);
                    }
                }
            }

            return slope;
        }
    }

    public static class DataSourceAlphaChi implements IDataSource, DataSourceIndependent {
        protected DataFunction chiData;
        protected DataDoubleArray alphaData;
        protected DataOverlap dsvo;
        protected final DataTag chiTag, alphaTag;
        protected DataInfoFunction chiInfo;
        protected DataInfoDoubleArray alphaInfo;
        protected final double errFac;
        protected int spiffiness;
        
        public DataSourceAlphaChi(DataOverlap dsvo, double errFac) {
            this.dsvo = dsvo;
            chiTag = new DataTag();
            alphaTag = new DataTag();
            this.errFac = errFac;
        }

        public void setSpiffiness(int newSpiffiness) {
            spiffiness = newSpiffiness;
        }

        public IData getData() {
            int numAlpha = dsvo.getAlphaSource().getNumAlpha();
            if (chiData == null || chiData.getLength() != numAlpha) {
                chiData = new DataFunction(new int[]{numAlpha});
                alphaData = new DataDoubleArray(chiData.getLength());
                chiInfo = new DataInfoFunction("chi", Null.DIMENSION, this);
                chiInfo.addTag(chiTag);
                alphaInfo = new DataInfoDoubleArray("alpha", Null.DIMENSION, new int[]{chiData.getLength()});
            }
            double[] chi = chiData.getData();
            for (int i=0; i<chi.length; i++) {
                if (spiffiness > 0 && dsvo.getRefSource() instanceof DataSourceOverlapAvgCollapsing) {
                    double avg = dsvo.getAverage(i);
                    double err = dsvo.getLogError(i);
                    chi[i] = avg * Math.exp(errFac*err);
                }
                else {
                    chi[i] = dsvo.getAverage(i) + errFac * dsvo.getError(i);
                    if (chi[i] < dsvo.getAverage(i)*0.01) {
                        chi[i] = dsvo.getAverage(i)*0.01;
                    }
                }
            }
            return chiData;
        }

        public DataTag getTag() {
            return chiTag;
        }

        public IDataInfo getDataInfo() {
            int numAlpha = dsvo.getAlphaSource().getNumAlpha();
            if (chiInfo == null || chiInfo.getLength() != numAlpha) {
                chiData = new DataFunction(new int[]{numAlpha});
                alphaData = new DataDoubleArray(chiData.getLength());
                chiInfo = new DataInfoFunction("chi", Null.DIMENSION, this);
                chiInfo.addTag(chiTag);
            }
            return chiInfo;
        }

        public DataDoubleArray getIndependentData(int j) {
            double[] alpha = alphaData.getData();
            for (int i=0; i<alpha.length; i++) {
                alpha[i] = dsvo.getAlphaSource().getAlpha(i);
            }
            return alphaData;
        }

        public DataInfoDoubleArray getIndependentDataInfo(int i) {
            int numAlpha = dsvo.getAlphaSource().getNumAlpha();
            if (alphaInfo == null || alphaInfo.getLength() != numAlpha) {
                alphaInfo = new DataInfoDoubleArray("alpha", Null.DIMENSION, new int[]{numAlpha});
                alphaInfo.addTag(alphaTag);
            }
            return alphaInfo;
        }

        public int getIndependentArrayDimension() {
            return 1;
        }

        public DataTag getIndependentTag() {
            return alphaTag;
        }
    }

    public static class DataSourceAlphaAlpha implements IDataSource, DataSourceIndependent {
        protected DataFunction chiData;
        protected DataDoubleArray alphaData;
        protected DataOverlap dsvo;
        protected final DataTag chiTag, alphaTag;
        protected DataInfoFunction chiInfo;
        protected DataInfoDoubleArray alphaInfo;
        
        public DataSourceAlphaAlpha(DataOverlap dsvo) {
            this.dsvo = dsvo;
            chiTag = new DataTag();
            alphaTag = new DataTag();
        }

        public IData getData() {
            int numAlpha = dsvo.getAlphaSource().getNumAlpha();
            if (chiData == null || chiData.getLength() != numAlpha) {
                chiData = new DataFunction(new int[]{numAlpha});
                alphaData = new DataDoubleArray(chiData.getLength());
                chiInfo = new DataInfoFunction("alpha", Null.DIMENSION, this);
                alphaInfo = new DataInfoDoubleArray("alpha", Null.DIMENSION, new int[]{chiData.getLength()});
            }
            double[] chi = chiData.getData();
            for (int i=0; i<chi.length; i++) {
                chi[i] = dsvo.getAlphaSource().getAlpha(i);
            }
            return chiData;
        }

        public DataTag getTag() {
            return chiTag;
        }

        public IDataInfo getDataInfo() {
            int numAlpha = dsvo.getAlphaSource().getNumAlpha();
            if (chiInfo == null || chiInfo.getLength() != numAlpha) {
                chiData = new DataFunction(new int[]{numAlpha});
                alphaData = new DataDoubleArray(chiData.getLength());
                chiInfo = new DataInfoFunction("alpha", Null.DIMENSION, this);
            }
            return chiInfo;
        }

        public DataDoubleArray getIndependentData(int j) {
            double[] alpha = alphaData.getData();
            for (int i=0; i<alpha.length; i++) {
                alpha[i] = dsvo.getAlphaSource().getAlpha(i);
            }
            return alphaData;
        }

        public DataInfoDoubleArray getIndependentDataInfo(int i) {
            int numAlpha = dsvo.getAlphaSource().getNumAlpha();
            if (alphaInfo == null || alphaInfo.getLength() != numAlpha) {
                alphaInfo = new DataInfoDoubleArray("alpha", Null.DIMENSION, new int[]{numAlpha});
            }
            return alphaInfo;
        }

        public int getIndependentArrayDimension() {
            return 1;
        }

        public DataTag getIndependentTag() {
            return alphaTag;
        }
    }
    
    public static class DataSourceUa extends DataSourceScalar {
        protected final DataOverlap dsvo;
        
        public DataSourceUa(DataOverlap dsvo) {
            super("Ua", Null.DIMENSION);
            this.dsvo = dsvo;
        }

        public double getDataAsScalar() {
            int numAlpha = dsvo.getAlphaSource().getNumAlpha();
            if (numAlpha < 4) {
                return Double.NaN;
            }

            double[] lnAlpha = new double[numAlpha];
            double[] lnAlphaDiff = new double[numAlpha];
            double[] lnRefAvg = new double[numAlpha];

            for (int j=0; j<numAlpha; j++) {
                double refOverlap = dsvo.getRefSource().getAverage(j);
                double targetOverlap = dsvo.getTargetSource().getAverage(j);
                lnRefAvg[j] = Math.log(targetOverlap);
                lnAlphaDiff[j] = Math.log(refOverlap/targetOverlap);

                double jAlpha = dsvo.getAlphaSource().getAlpha(j);
                lnAlpha[j] = Math.log(jAlpha);
                lnAlphaDiff[j] -= lnAlpha[j];
            }

            double lnUa = Double.NaN;
            if (lnAlphaDiff[0] < 0) {
                // first new alpha is less than initial first alpha
                return Double.NaN;
            }
            else if (lnAlphaDiff[numAlpha-1] > 0) {
                return Double.NaN;
            }
            AkimaSpline spline = new AkimaSpline();

            double[] x = new double[]{Math.log(dsvo.getOverlapAverageAndError()[0])};
            spline.setInputData(lnAlpha, lnRefAvg);
            lnUa = spline.doInterpolation(x)[0];

            return 2*Math.exp(lnUa+x[0]);
        }
    }

    public static class DataSourceUa2 extends DataSourceScalar {

        protected final IntegratorOverlap integratorOS;
        protected final DataOverlap dsvo;
        protected final AccumulatorAverage accRefUa2, accTargetUa2;
        
        public DataSourceUa2(IntegratorOverlap integratorOS, DataOverlap dsvo, AccumulatorAverage accRef, AccumulatorAverage accTarget) {
            super("Ua2", Null.DIMENSION);
            this.integratorOS = integratorOS;
            this.dsvo = dsvo;
            this.accRefUa2 = accRef;
            this.accTargetUa2 = accTarget;
        }

        public double getDataAsScalar() {
            int numAlpha = dsvo.getAlphaSource().getNumAlpha();
            if (numAlpha == 1) {
                return Double.NaN;
            }
            

            double[] lnAlpha = new double[numAlpha];
            double[] lnRefAvg = new double[numAlpha];
            double[] lnTargetAvg = new double[numAlpha];

            IData refData = ((DataGroup)accRefUa2.getData()).getData(AccumulatorAverage.AVERAGE.index);
            IData targetData = ((DataGroup)accTargetUa2.getData()).getData(AccumulatorAverage.AVERAGE.index);
            
            for (int j=0; j<numAlpha; j++) {
                double jAlpha = dsvo.getAlphaSource().getAlpha(j);

                lnRefAvg[j] = Math.log(refData.getValue(j));
                lnTargetAvg[j] = Math.log(targetData.getValue(j));

                lnAlpha[j] = Math.log(jAlpha);
            }

            double[] x = new double[1];
            x[0] = Math.log(dsvo.getOverlapAverageAndError()[0]);
            if (x[0] <= lnAlpha[0]) {
                // first new alpha is less than initial first alpha
                return Double.NaN;
            }
            else if (x[0] >= lnAlpha[numAlpha-1]) {
                return Double.NaN;
            }
            double refUa2 = Double.NaN;
            double targetUa2 = Double.NaN;
            if (numAlpha > 4) {
                AkimaSpline spline = new AkimaSpline();
                spline.setInputData(lnAlpha, lnRefAvg);
                refUa2 = spline.doInterpolation(x)[0];
                spline.setInputData(lnAlpha, lnTargetAvg);
                targetUa2 = spline.doInterpolation(x)[0];
            }
            else {
                //linear interpolation (only 3 points)
                double myLnAlpha = x[0];
                for (int i=0; i<numAlpha-1; i++) {
                    if (lnAlpha[i+1] > myLnAlpha) {
                        refUa2 = lnRefAvg[i] + (lnRefAvg[i+1] - lnRefAvg[i]) / (lnAlpha[i+1]-lnAlpha[i])*(myLnAlpha-lnAlpha[i]);
                        targetUa2 = lnTargetAvg[i] + (lnTargetAvg[i+1] - lnTargetAvg[i]) / (lnAlpha[i+1]-lnAlpha[i])*(myLnAlpha-lnAlpha[i]);
                    }
                }
            }

            return 2*Math.exp(refUa2) + 2*Math.exp(2*x[0]+targetUa2);
        }
    }

    public static class DataSinkReweightedDeltaU implements IDataSink {
        protected double[] sum, sumWeights, ratio;
        protected double[] lnAlpha;
        protected final AlphaSource alphaSource;
        protected final AkimaSpline akima;
        
        public DataSinkReweightedDeltaU(AlphaSource alphaSource) {
            this.alphaSource = alphaSource;
            akima = new AkimaSpline();
            lnAlpha = new double[0];
            ratio = new double[0];
        }
        
        public void reset() {
            int numAlpha = alphaSource.getNumAlpha();
            sum = new double[numAlpha];
            sumWeights = new double[numAlpha];
        }

        public double getDeltaU(double alpha) {
            if (lnAlpha.length != alphaSource.getNumAlpha()) {
                lnAlpha = new double[alphaSource.getNumAlpha()];
                ratio = new double[alphaSource.getNumAlpha()];
            }
            for (int i=0; i<lnAlpha.length; i++) {
                lnAlpha[i] = Math.log(alphaSource.getAlpha(i));
                ratio[i] = sum[i]/sumWeights[i];
            }
            if (lnAlpha.length == 1) return ratio[0];
            akima.setInputData(lnAlpha, ratio);
            return akima.doInterpolation(new double[]{Math.log(alpha)})[0];
        }

        public void putData(IData data) {
            for (int i=0; i<data.getLength(); i++) {
                double x = data.getValue(i);
                if (x > 0) {
                    sum[i] += x*Math.log(x);
                    sumWeights[i] += x;
                }
            }
        }

        public void putDataInfo(IDataInfo dataInfo) {
            sum = new double[dataInfo.getLength()];
            sumWeights = new double[sum.length];
        }
    }

    public static class Applet extends javax.swing.JApplet {

        public void init() {
            final MultiharmonicMC sim = new MultiharmonicMC();
            MultiharmonicGraphicMC simGraphic = new MultiharmonicGraphicMC(sim, sim.getSpace());
            getContentPane().add(simGraphic.getPanel());
        }

        private static final long serialVersionUID = 1L;
    }

}
