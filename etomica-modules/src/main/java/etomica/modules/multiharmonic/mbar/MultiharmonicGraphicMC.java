/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic.mbar;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate2;
import etomica.atom.DiameterHashByType;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.graphics.*;
import etomica.math.function.Function;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.space.Space;
import etomica.space1d.Vector1D;
import etomica.units.Pixel;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;
import net.miginfocom.swing.MigLayout;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;

public class MultiharmonicGraphicMC extends SimulationGraphic {

    private final static String APP_NAME = "MultiharmonicMC";
    protected final static int REPAINT_INTERVAL = 300000;

    protected final MultiharmonicMC sim;

    /**
     * 
     */
    public MultiharmonicGraphicMC(MultiharmonicMC simulation, Space _space) {
        super(simulation, GRAPHIC_ONLY, APP_NAME, REPAINT_INTERVAL);
        this.sim = simulation;
        final DisplayBox displayBoxA = getDisplayBox(sim.boxA);
        remove(displayBoxA);
        final DisplayBox displayBoxB = getDisplayBox(sim.boxB);
        remove(displayBoxB);
        final Device controllerButtons = getController();
        getPanel().graphicsPanel.remove(controllerButtons.graphic());
        getPanel().footerPanel.add(controllerButtons.graphic());
        getPanel().graphicsPanel.setLayout(new MigLayout());

        displayBoxA.setPixelUnit(new Pixel(350/sim.boxA.getBoundary().getBoxSize().getX(0)));
        ((DiameterHashByType)displayBoxA.getDiameterHash()).setDiameter(sim.species.getLeafType(), 0.02);
        displayBoxB.setPixelUnit(new Pixel(350/sim.boxB.getBoundary().getBoxSize().getX(0)));
        displayBoxB.setDiameterHash(displayBoxA.getDiameterHash());

        final DisplayPlotXChart plot = new DisplayPlotXChart();
        
        DataSourceCountSteps stepCounter = new DataSourceCountSteps(sim.integratorOS);

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

        final boolean isBroken = false;
        
        DataSourceScalar myDataSourceMin = new DataSourceFE("average", sim.meterOverlapA, sim.meterOverlapB, isBroken ? -1 : 0);
        AccumulatorHistory historyMin = new AccumulatorHistory(new HistoryCollapsingDiscard(102, 3));
        DataPumpListener myPump = new DataPumpListener(myDataSourceMin, historyMin, 10);
        sim.integratorOS.getEventManager().addListener(myPump);
        historyMin.setTimeDataSource(stepCounter);
        historyMin.setPushInterval(10);
        historyMin.setDataSink(plot.getDataSet().makeDataSink());
        plot.setLegend(new DataTag[]{historyMin.getTag()}, isBroken ? "min" : "measured");
        if (isBroken) {
            DataSourceScalar myDataSourceMax = new DataSourceFE("average", sim.meterOverlapA, sim.meterOverlapB, 1);
            AccumulatorHistory historyMax = new AccumulatorHistory(new HistoryCollapsingDiscard(102, 3));
            myPump = new DataPumpListener(myDataSourceMax, historyMax, 10);
            sim.integratorOS.getEventManager().addListener(myPump);
            historyMax.setTimeDataSource(stepCounter);
            historyMax.setPushInterval(10);
            historyMax.setDataSink(plot.getDataSet().makeDataSink());
            plot.setLegend(new DataTag[]{historyMax.getTag()}, "max");
        }

        DisplayPlotXChart alphaChiPlot = new DisplayPlotXChart();
        alphaChiPlot.getPlot().setXLog(true);
        alphaChiPlot.getPlot().setYLog(true);
        IDataSource alphaChi = new DataSourceAlphaChi(sim.meterOverlapA, sim.meterOverlapB, isBroken);
        DataPumpListener alphaChiPump = new DataPumpListener(alphaChi, alphaChiPlot.getDataSet().makeDataSink(), 10);
        sim.integratorOS.getEventManager().addListener(alphaChiPump);
        alphaChiPlot.setDoLegend(false);
        alphaChiPlot.getPlot().setTitle("chi / alpha");
        
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
        x0Slider.setValue(0.0);
        omegaASlider.setModifier(new ModifierGeneral(sim.potentialA, "springConstant"));
        omegaASlider.setMinimum(0.1);
        omegaASlider.setMaximum(50.0);
        omegaASlider.setValue(1.0);
        omegaBSlider.setModifier(new ModifierGeneral(sim.potentialB, "springConstant"));
        omegaBSlider.setMinimum(0.1);
        omegaBSlider.setMaximum(10.0);
        omegaBSlider.setValue(1.0);
        
        DeviceBox alphaCenterBox = new DeviceBox();
        alphaCenterBox.setController(sim.getController());
        alphaCenterBox.setModifier(new Modifier() {
            
            public void setValue(double newValue) {
                if (newValue <= 0) throw new IllegalArgumentException();
                double[] span = sim.meterOverlapA.getAlphaSpan();
                sim.meterOverlapA.setAlpha(new double[]{newValue}, span);
                sim.meterOverlapB.setAlpha(new double[]{newValue}, span);
            }
            
            public double getValue() {
                return sim.meterOverlapA.getAlphaCenter()[0];
            }
            
            public String getLabel() {
                return "alpha center";
            }
            
            public Dimension getDimension() {
                return Null.DIMENSION;
            }
        });
        DeviceBox alphaSpanBox = new DeviceBox();
        alphaSpanBox.setController(sim.getController());
        alphaSpanBox.setModifier(new Modifier() {
            
            public void setValue(double newValue) {
                if (newValue <= 0 || newValue >= 50) throw new IllegalArgumentException();
                double[] center = sim.meterOverlapA.getAlphaCenter();
                sim.meterOverlapA.setAlpha(center, new double[]{newValue});
                sim.meterOverlapB.setAlpha(center, new double[]{newValue});
            }
            
            public double getValue() {
                return sim.meterOverlapA.getAlphaSpan()[0];
            }
            
            public String getLabel() {
                return "alpha span";
            }
            
            public Dimension getDimension() {
                return Null.DIMENSION;
            }
        });
        
        DeviceBox numAlphaBox = new DeviceBox();
        numAlphaBox.setController(sim.getController());
        numAlphaBox.setInteger(true);
        numAlphaBox.setModifier(new Modifier() {
            
            public void setValue(double newValue) {
                if (newValue <= 0 || newValue >= 50 || newValue%2 == 0) throw new IllegalArgumentException();
                sim.meterOverlapA.setNumAlpha((int)newValue);
                sim.meterOverlapB.setNumAlpha((int)newValue);
            }
            
            public double getValue() {
                return sim.meterOverlapA.getNumAlpha();
            }
            
            public String getLabel() {
                return "# alpha";
            }
            
            public Dimension getDimension() {
                return Null.DIMENSION;
            }
        });

        
        DeviceNSelector nSlider = new DeviceNSelector(sim.getController());
        nSlider.setBox(sim.boxA);
        nSlider.setSpecies(sim.species);
        nSlider.setMaximum(50);
        nSlider.setMinimum(1);
        nSlider.setShowValues(true);
        nSlider.setLabel("Number of atoms");
        nSlider.setShowBorder(true);

        DataSourceScalar delta = new DataSourceScalar("exact",Energy.DIMENSION) {
            public double getDataAsScalar() {
                return 0.5*sim.boxA.getLeafList().size() * Math.log(omegaBSlider.getValue()/omegaASlider.getValue());
            }
        };
        
        AccumulatorHistory deltaHistory = new AccumulatorHistory(new HistoryCollapsingDiscard(102, 3));
        deltaHistory.setPushInterval(10);
        DataPumpListener exactPump = new DataPumpListener(delta, deltaHistory, 10);
        deltaHistory.setDataSink(plot.getDataSet().makeDataSink());
        sim.integratorOS.getEventManager().addListener(exactPump);
        dataStreamPumps.add(exactPump);
        deltaHistory.setTimeDataSource(stepCounter);
        
        final DisplayPlotXChart uPlot = new DisplayPlotXChart();
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
        uA.getXSource().setXMin(-sim.boxA.getBoundary().getBoxSize().getX(0) / 2);
        uB.getXSource().setXMin(-sim.boxA.getBoundary().getBoxSize().getX(0) / 2);
        uA.getXSource().setXMax(sim.boxA.getBoundary().getBoxSize().getX(0) / 2);
        uB.getXSource().setXMax(sim.boxA.getBoundary().getBoxSize().getX(0) / 2);
        final DataPump uAPump = new DataPump(uA, uPlot.getDataSet().makeDataSink());
        final DataPump uBPump = new DataPump(uB, uPlot.getDataSet().makeDataSink());
        IAction uUpdate = new IAction() {
            public void actionPerformed() {
                uA.update();
                uB.update();
                uAPump.actionPerformed();
                uBPump.actionPerformed();
            }
        };
        omegaASlider.setPostAction(uUpdate);
        omegaBSlider.setPostAction(uUpdate);
        x0Slider.setPostAction(uUpdate);

        uPlot.getDataSet().setUpdatingOnAnyChange(true);

        //controls -- start/pause and sliders
        JTabbedPane sliderPanel = new JTabbedPane();
        sliderPanel.add(x0Slider.graphic(), "x0");
        sliderPanel.add(omegaASlider.graphic(), "omegaA");
        sliderPanel.add(omegaBSlider.graphic(), "omegaB");
        
        JPanel topPanel = new JPanel(new MigLayout());
        JPanel displayPanel = new JPanel(new MigLayout("flowy, novisualpadding"));
        JPanel controlsPanel = new JPanel(new MigLayout("flowy, novisualpadding"));

        controlsPanel.add(sliderPanel, "align center");
        controlsPanel.add(nSlider.graphic(), "align center");

        JPanel alphaPanel = new JPanel(new MigLayout());
        alphaPanel.add(alphaCenterBox.graphic());
        alphaPanel.add(alphaSpanBox.graphic());
        alphaPanel.add(numAlphaBox.graphic());
        controlsPanel.add(alphaPanel, "align center");

        displayPanel.add(uPlot.graphic());

        JTabbedPane displayBoxPanel = new JTabbedPane();
        displayBoxPanel.add(displayBoxA.graphic(), "box A");
        displayBoxPanel.add(displayBoxB.graphic(), "box B");
        displayPanel.add(displayBoxPanel, "growx, center, gapx 4% 9%");
        topPanel.add(controlsPanel, "");
        topPanel.add(displayPanel, "");
        
        getPanel().graphicsPanel.setLayout(new MigLayout());
        getPanel().graphicsPanel.add(topPanel, "wrap");

        JTabbedPane plotTabs = new JTabbedPane();

        getPanel().graphicsPanel.add(plotTabs);

        getController().getReinitButton().setPostAction(new IAction() {
        	public void actionPerformed() {
                displayBoxA.repaint();
                displayBoxB.repaint();
                plot.getPanel().repaint();
        	}
        });

        uUpdate.actionPerformed();

//        AccumulatorAverageCollapsingLog accumulatorEnergy = new AccumulatorAverageCollapsingLog(sim.getRandom());
//        DataPump dataPumpEnergy = new DataPump(sim.meter, accumulatorEnergy);
//        dataStreamPumps.add(dataPumpEnergy);
//
//        DataFork feFork = new DataFork();
//        DataProcessorFunction negative = new DataProcessorFunction(new Function() {
//            public double f(double x) {
//                return -x;
//            }
//        });
//        accumulatorEnergy.setDataSink(negative);
//        negative.setDataSink(feFork);
//        DisplayPlot funPlot = new DisplayPlot();
//        funPlot.getPlot().setTitle("Free Energy Difference Convergence");
//        funPlot.setDoLegend(false);
//        funPlot.setSize(350, 250);
//        feFork.addDataSink(funPlot.getDataSet().makeDataSink());
//        funPlot.getPlot().setXLog(true);
//        
//        DataProcessorBounds dpUpper = new DataProcessorBounds();
//        dpUpper.setAccumulator(accumulatorEnergy);
//        dpUpper.setIsUpperBound(true);
//        feFork.addDataSink(dpUpper);
//        dpUpper.setDataSink(funPlot.getDataSet().makeDataSink());
//        funPlot.setLegend(new DataTag[]{dpUpper.getTag()}, "FE upper bound");
//        DataProcessorBounds dpLower = new DataProcessorBounds();
//        dpLower.setAccumulator(accumulatorEnergy);
//        dpLower.setIsUpperBound(false);
//        feFork.addDataSink(dpLower);
//        dpLower.setDataSink(funPlot.getDataSet().makeDataSink());
//        funPlot.setLegend(new DataTag[]{dpLower.getTag()}, "FE lower bound");
//        
//        tab2.add(funPlot.graphic());
        plotTabs.add(plot.graphic(), "Free Energy");
        plotTabs.add(alphaChiPlot.graphic(), "Chi");
//
//        DataProcessorVar dpstdev = new DataProcessorVar();
//        dpstdev.setAccumulator(accumulatorEnergy);
//        feFork.addDataSink(dpstdev);
//        DataFork stdevFork = new DataFork();
//        dpstdev.setDataSink(stdevFork);
//        DataProcessorFunction lnstdev = new DataProcessorFunction(new Function.Log());
//        stdevFork.addDataSink(lnstdev);
//
//        DataProcessorDyDLnx dAdlnx = new DataProcessorDyDLnx();
//        feFork.addDataSink(dAdlnx);
//        
//        DataProcessorFunction ndAdlnx = new DataProcessorFunction(new Function() {
//            public double f(double x) {return -x;}
//        });
//        dAdlnx.setDataSink(ndAdlnx);
//        DataFork ndAdlnxFork = new DataFork();
//        ndAdlnx.setDataSink(ndAdlnxFork);
//        DataProcessorFunction lndAdlnx = new DataProcessorFunction(new Function.Log());
//        ndAdlnxFork.addDataSink(lndAdlnx);
//        DataProcessorDyDLnx dlndAdlnx = new DataProcessorDyDLnx();
//        lndAdlnx.setDataSink(dlndAdlnx);
//
//        DisplayPlot dfunPlot = new DisplayPlot();
//        dfunPlot.getPlot().setTitle("d(ln(x))/dlnN");
//        dfunPlot.setSize(350, 250);
//        dlndAdlnx.setDataSink(dfunPlot.getDataSet().makeDataSink());
//        dfunPlot.setLegend(new DataTag[]{dlndAdlnx.getTag()}, "x=-dA/dlnN");
//        dfunPlot.getPlot().setXLog(true);
//        dfunPlot.setDoLegend(true);
//
//        DataProcessorBias dpBias = new DataProcessorBias(delta);
//        feFork.addDataSink(dpBias);
//        DataFork biasFork = new DataFork();
//        dpBias.setDataSink(biasFork);
//        DataProcessorFunction lnb = new DataProcessorFunction(new Function.Log());
//        biasFork.addDataSink(lnb);
//        
//        DataProcessorDyDLnx dlnbDlnx = new DataProcessorDyDLnx();
//        lnb.setDataSink(dlnbDlnx);
//        
//        dlnbDlnx.setDataSink(dfunPlot.getDataSet().makeDataSink());
//        dfunPlot.setLegend(new DataTag[]{dlnbDlnx.getTag()}, "x=b");
//        
//        DataProcessorDyDLnx dlnstdevDlnx = new DataProcessorDyDLnx();
//        lnstdev.setDataSink(dlnstdevDlnx);
//        dlnstdevDlnx.setDataSink(dfunPlot.getDataSet().makeDataSink());
//        dfunPlot.setLegend(new DataTag[]{dlnstdevDlnx.getTag()}, "x=var_A");
//
//        DisplayPlot lfunPlot = new DisplayPlot();
//        lfunPlot.setSize(300, 250);
//        lfunPlot.getPlot().setTitle("x");
//        ndAdlnxFork.addDataSink(lfunPlot.getDataSet().makeDataSink());
//        lfunPlot.setLegend(new DataTag[]{ndAdlnxFork.getTag()}, "x=-dA/dlnN");
//        lfunPlot.getPlot().setXLog(true);
//        lfunPlot.getPlot().setYLog(true);
//
//        stdevFork.addDataSink(lfunPlot.getDataSet().makeDataSink());
//        lfunPlot.setLegend(new DataTag[]{stdevFork.getTag()}, "x=var_A");
//        
//        biasFork.addDataSink(lfunPlot.getDataSet().makeDataSink());
//        lfunPlot.setLegend(new DataTag[]{biasFork.getTag()}, "x=b");
//
//        JPanel tab3 = new JPanel(new GridBagLayout());
//        plotTabs.add(tab3, "Bias");
//        tab3.add(lfunPlot.graphic());
//        tab3.add(dfunPlot.graphic(), horizGBC);

        sim.getController2().addActivity(new ActivityIntegrate2(sim.integratorOS)).setSleepPeriod(0);
    }

    public static void main(String[] args) {
        final MultiharmonicMC sim = new MultiharmonicMC();
        MultiharmonicGraphicMC simGraphic = new MultiharmonicGraphicMC(sim, sim.getSpace());
        SimulationGraphic.makeAndDisplayFrame(simGraphic.getPanel(), APP_NAME);
    }

    public static class DataSourceFE extends DataSourceScalar {
        protected final MeterMBAR meterA, meterB;
        protected final int broken;

        public DataSourceFE(String label, MeterMBAR meterA, MeterMBAR meterB, int broken) {
            super(label, Null.DIMENSION);
            this.meterA = meterA;
            this.meterB = meterB;
            this.broken = broken;;
        }

        public double getDataAsScalar() {
            DataDoubleArray dataA = (DataDoubleArray)meterA.getData();
            DataDoubleArray dataB = (DataDoubleArray)meterB.getData();
            double prevDelta = 0;
            double prevAlpha = 0;
            int n = meterA.getNumAlpha();
            int i = 0;
            if (broken == 1) {
                // iterate backwards
                i = n-1;
            }
                
            while (true) {
                double sumAA = dataA.getValue(new int[]{i,0});
                double sumAB = dataA.getValue(new int[]{i,1});
                double sumBA = dataB.getValue(new int[]{i,0});
                double sumBB = dataB.getValue(new int[]{i,1});
                double lnAlpha = Math.log(meterA.getAlpha(i)[0]);
                double delta;
                if (broken == 0) {
                    delta = Math.log(-sumBB/sumAB);
                }
                else {
                    sumAA += meterA.getCallCount();
                    sumBB += meterB.getCallCount() / meterB.getAlpha(i)[0];
                    double lnchi = -Math.log((sumAB + sumBB) / (sumAA + sumBA));
                    delta = lnchi - lnAlpha;
                }
                if (((i>0 && broken != 1) || (i<n-1 && broken==1)) && prevDelta*delta<=0) {
                    return prevAlpha +  (lnAlpha - prevAlpha) / (delta - prevDelta) * (-prevDelta);
                }
                prevDelta = delta;
                prevAlpha = lnAlpha;
                if (broken == 1) {
                    i--;
                    if (i==-1) break;
                }
                else {
                    i++;
                    if (i==n) break;
                }
            }
            return Double.NaN;
        }
    }

    public static class DataSourceAlphaChi implements IDataSource, DataSourceIndependent {
        protected DataFunction chiData;
        protected DataDoubleArray alphaData;
        protected MeterMBAR meterOverlapA, meterOverlapB;
        protected final DataTag chiTag, alphaTag;
        protected DataInfoFunction chiInfo;
        protected DataInfoDoubleArray alphaInfo;
        protected boolean isBroken;
        
        public DataSourceAlphaChi(MeterMBAR meterOverlapA, MeterMBAR meterOverlapB, boolean isBroken) {
            this.meterOverlapA = meterOverlapA;
            this.meterOverlapB = meterOverlapB;
            chiTag = new DataTag();
            alphaTag = new DataTag();
            this.isBroken = isBroken;
        }

        public IData getData() {
            if (chiData == null || chiData.getLength() != meterOverlapA.getNumAlpha()) {
                chiData = new DataFunction(new int[]{meterOverlapA.getNumAlpha()});
                alphaData = new DataDoubleArray(chiData.getLength());
                chiInfo = new DataInfoFunction("chi", Null.DIMENSION, this);
                alphaInfo = new DataInfoDoubleArray("alpha", Null.DIMENSION, new int[]{chiData.getLength()});
            }
            double[] chi = chiData.getData();
            for (int i=0; i<chi.length; i++) {
                DataDoubleArray dataA = (DataDoubleArray)meterOverlapA.getData();
                DataDoubleArray dataB = (DataDoubleArray)meterOverlapB.getData();
                double sumAA = dataA.getValue(new int[]{i,0});
                double sumAB = dataA.getValue(new int[]{i,1});
                double sumBA = dataB.getValue(new int[]{i,0});
                double sumBB = dataB.getValue(new int[]{i,1});
//                System.out.println(i+" "+sumAA+" "+sumAB+" "+sumBA+" "+sumBB);
                if (isBroken) {
                    chi[i] = (sumAB + sumBB) / (sumAA + sumBA) / meterOverlapA.getAlpha(i)[1];
                    sumAA += meterOverlapA.getCallCount();
                    sumBB += meterOverlapB.getCallCount() / meterOverlapB.getAlpha(i)[0];
                }
                else {
                    chi[i] = meterOverlapA.getAlpha(i)[0] / ((-sumAB/meterOverlapA.getCallCount())/(sumBB/meterOverlapB.getCallCount()));
                }
            }
            return chiData;
        }

        public DataTag getTag() {
            return chiTag;
        }

        public IDataInfo getDataInfo() {
            if (chiInfo == null || chiInfo.getLength() != meterOverlapA.getNumAlpha()) {
                chiData = new DataFunction(new int[]{meterOverlapA.getNumAlpha()});
                alphaData = new DataDoubleArray(chiData.getLength());
                chiInfo = new DataInfoFunction("chi", Null.DIMENSION, this);
            }
            return chiInfo;
        }

        public DataDoubleArray getIndependentData(int j) {
            double[] alpha = alphaData.getData();
            for (int i=0; i<alpha.length; i++) {
                alpha[i] = meterOverlapA.getAlpha(i)[0];
            }
            return alphaData;
        }

        public DataInfoDoubleArray getIndependentDataInfo(int i) {
            if (alphaInfo == null || alphaInfo.getLength() != meterOverlapA.getNumAlpha()) {
                alphaInfo = new DataInfoDoubleArray("alpha", Null.DIMENSION, new int[]{meterOverlapA.getNumAlpha()});
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


    public static class Applet extends javax.swing.JApplet {

        public void init() {
            final MultiharmonicMC sim = new MultiharmonicMC();
            MultiharmonicGraphicMC simGraphic = new MultiharmonicGraphicMC(sim, sim.getSpace());
            getContentPane().add(simGraphic.getPanel());
        }

        private static final long serialVersionUID = 1L;
    }

}
