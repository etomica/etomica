/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic.umbrella;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate2;
import etomica.atom.DiameterHashByType;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.data.types.DataGroup;
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
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;
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
        final DisplayBox displayBoxA = getDisplayBox(sim.box);
        remove(displayBoxA);
        final Device controllerButtons = getController();
        getPanel().graphicsPanel.remove(controllerButtons.graphic());
        getPanel().footerPanel.add(controllerButtons.graphic());
        getPanel().graphicsPanel.setLayout(new MigLayout("flowy"));
        
        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();
        dataStreamPumps.add(simulation.dataPumpA);

        ((DiameterHashByType)displayBoxA.getDiameterHash()).setDiameter(sim.species.getLeafType(), 0.02);

        final DisplayPlotXChart plot = new DisplayPlotXChart();
        
        DataSourceCountSteps stepCounter = new DataSourceCountSteps(sim.integrator);

        DataProcessorFunction log = new DataProcessorFunction(new Function() {
            public double f(double x) {return -Math.log(x);}
        });
        DataSourceScalar myDataSource = new DataSourceScalar("average", Null.DIMENSION) {
            public double getDataAsScalar() {
                return ((DataGroup)sim.accumulator.getData()).getData(sim.accumulator.RATIO.index).getValue(1);
            }
        };
        DataPumpListener myPump = new DataPumpListener(myDataSource, log, 10);
        sim.integrator.getEventManager().addListener(myPump);
        AccumulatorHistory history = new AccumulatorHistory(new HistoryCollapsingDiscard(102, 3));
        history.setTimeDataSource(stepCounter);
        log.setDataSink(history);
        history.setPushInterval(10);
        history.setDataSink(plot.getDataSet().makeDataSink());
        
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
        
        DeviceBox alphaFacBox = new DeviceBox();
        alphaFacBox.setModifier(new Modifier() {
            
            public void setValue(double newValue) {
                if (newValue <= 0) throw new IllegalArgumentException();
                sim.meterUmbrella.setAlphaFac(newValue);
                sim.integrator.getMoveManager().setFrequency(sim.moveA, newValue);
                sim.accumulator.reset();
            }
            
            public double getValue() {
                return sim.meterUmbrella.getAlphaFac();
            }
            
            public String getLabel() {
                return "alpha fac";
            }
            
            public Dimension getDimension() {
                return Null.DIMENSION;
            }
        });
        alphaFacBox.doUpdate();
        
        DeviceNSelector nSlider = new DeviceNSelector(sim.getController());
        nSlider.setBox(sim.box);
        nSlider.setSpecies(sim.species);
        nSlider.setMaximum(50);
        nSlider.setMinimum(1);
        nSlider.setShowValues(true);
        nSlider.setLabel("Number of atoms");
        nSlider.setShowBorder(true);

        DataSourceScalar delta = new DataSourceScalar("exact",Energy.DIMENSION) {
            public double getDataAsScalar() {
                return 0.5*sim.box.getLeafList().size() * Math.log(omegaBSlider.getValue()/omegaASlider.getValue());
            }
        };
        
        AccumulatorHistory deltaHistory = new AccumulatorHistory(new HistoryCollapsingDiscard(102, 3));
        deltaHistory.setPushInterval(10);
        DataPumpListener exactPump = new DataPumpListener(delta, deltaHistory, 10);
        deltaHistory.setDataSink(plot.getDataSet().makeDataSink());
        sim.integrator.getEventManager().addListener(exactPump);
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
        uA.getXSource().setXMin(-sim.box.getBoundary().getBoxSize().getX(0) / 2);
        uB.getXSource().setXMin(-sim.box.getBoundary().getBoxSize().getX(0) / 2);
        uA.getXSource().setXMax(sim.box.getBoundary().getBoxSize().getX(0) / 2);
        uB.getXSource().setXMax(sim.box.getBoundary().getBoxSize().getX(0) / 2);
        final DataPump uAPump = new DataPump(uA, uPlot.getDataSet().makeDataSink());
        final DataPump uBPump = new DataPump(uB, uPlot.getDataSet().makeDataSink());
        IAction uUpdate = new IAction() {
            public void actionPerformed() {
                uA.update();
                uB.update();
                uAPump.actionPerformed();
                uBPump.actionPerformed();
                double alpha = Math.pow(omegaBSlider.getValue()/omegaASlider.getValue(), -0.5*sim.box.getLeafList().size());
                sim.meterUmbrella.setAlpha(alpha);
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
        
        JPanel topPanel = new JPanel(new MigLayout("flowy, novisualpadding"));
        JPanel displayPanel = new JPanel(new MigLayout("flowy"));
        topPanel.add(sliderPanel, "");
        topPanel.add(nSlider.graphic(), "");
        topPanel.add(alphaFacBox.graphic(), "aligny top");
        displayPanel.add(uPlot.graphic(), "");
        displayPanel.add(displayBoxA.graphic(), "growx, center, gapx 4% 9%");
        topPanel.add(displayPanel, "newline, spany");
        
        getPanel().graphicsPanel.add(topPanel);

        JTabbedPane plotTabs = new JTabbedPane();

        getPanel().graphicsPanel.add(plotTabs);

        getController().getReinitButton().setPostAction(new IAction() {
        	public void actionPerformed() {
                displayBoxA.repaint();
                plot.getPlot().doUpdate();
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
//        DisplayPlotXChart funPlot = new DisplayPlotXChart();
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
        JPanel tab2 = new JPanel(new MigLayout());
        plotTabs.add(tab2, "Free Energy");
//        tab2.add(funPlot.graphic());
        tab2.add(plot.graphic());
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
//        DisplayPlotXChart dfunPlot = new DisplayPlotXChart();
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
//        DisplayPlotXChart lfunPlot = new DisplayPlotXChart();
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
//        JPanel tab3 = new JPanel(new MigLayout());
//        plotTabs.add(tab3, "Bias");
//        tab3.add(lfunPlot.graphic());
//        tab3.add(dfunPlot.graphic(), horizGBC);

        sim.getController2().addActivity(new ActivityIntegrate2(sim.integrator)).setSleepPeriod(0);
    }

    public static void main(String[] args) {
        final MultiharmonicMC sim = new MultiharmonicMC();
        MultiharmonicGraphicMC simGraphic = new MultiharmonicGraphicMC(sim, sim.getSpace());
        SimulationGraphic.makeAndDisplayFrame(simGraphic.getPanel(), APP_NAME);
    }

    public static class DataSourceAlphaChi implements IDataSource, DataSourceIndependent {
        protected DataFunction chiData;
        protected DataDoubleArray alphaData;
        protected DataSourceVirialOverlap dsvo;
        protected final DataTag chiTag, alphaTag;
        protected DataInfoFunction chiInfo;
        protected DataInfoDoubleArray alphaInfo;
        protected final double errFac;
        
        public DataSourceAlphaChi(DataSourceVirialOverlap dsvo, double errFac) {
            this.dsvo = dsvo;
            chiTag = new DataTag();
            alphaTag = new DataTag();
            this.errFac = errFac;
        }

        public IData getData() {
            if (chiData == null || chiData.getLength() != dsvo.getAccumulators()[0].getNBennetPoints()) {
                chiData = new DataFunction(new int[]{dsvo.getAccumulators()[0].getNBennetPoints()});
                alphaData = new DataDoubleArray(chiData.getLength());
                chiInfo = new DataInfoFunction("chi", Null.DIMENSION, this);
                alphaInfo = new DataInfoDoubleArray("alpha", Null.DIMENSION, new int[]{chiData.getLength()});
            }
            double[] chi = chiData.getData();
            for (int i=0; i<chi.length; i++) {
                chi[i] = dsvo.getAverage(i) + errFac * dsvo.getError(i);
                if (chi[i] < dsvo.getAverage(i)*0.01) {
                    chi[i] = dsvo.getAverage(i)*0.01;
                }
            }                
            return chiData;
        }

        public DataTag getTag() {
            return chiTag;
        }

        public IDataInfo getDataInfo() {
            if (chiInfo == null || chiInfo.getLength() != dsvo.getAccumulators()[0].getNBennetPoints()) {
                chiData = new DataFunction(new int[]{dsvo.getAccumulators()[0].getNBennetPoints()});
                alphaData = new DataDoubleArray(chiData.getLength());
                chiInfo = new DataInfoFunction("chi", Null.DIMENSION, this);
            }
            return chiInfo;
        }

        public DataDoubleArray getIndependentData(int j) {
            double[] alpha = alphaData.getData();
            AccumulatorVirialOverlapSingleAverage acc = dsvo.getAccumulators()[0];
            for (int i=0; i<alpha.length; i++) {
                alpha[i] = acc.getBennetBias(i);
            }
            return alphaData;
        }

        public DataInfoDoubleArray getIndependentDataInfo(int i) {
            if (alphaInfo == null || alphaInfo.getLength() != dsvo.getAccumulators()[0].getNBennetPoints()) {
                alphaInfo = new DataInfoDoubleArray("alpha", Null.DIMENSION, new int[]{dsvo.getAccumulators()[0].getNBennetPoints()});
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
        protected DataSourceVirialOverlap dsvo;
        protected final DataTag chiTag, alphaTag;
        protected DataInfoFunction chiInfo;
        protected DataInfoDoubleArray alphaInfo;
        
        public DataSourceAlphaAlpha(DataSourceVirialOverlap dsvo) {
            this.dsvo = dsvo;
            chiTag = new DataTag();
            alphaTag = new DataTag();
        }

        public IData getData() {
            if (chiData == null || chiData.getLength() != dsvo.getAccumulators()[0].getNBennetPoints()) {
                chiData = new DataFunction(new int[]{dsvo.getAccumulators()[0].getNBennetPoints()});
                alphaData = new DataDoubleArray(chiData.getLength());
                chiInfo = new DataInfoFunction("alpha", Null.DIMENSION, this);
                alphaInfo = new DataInfoDoubleArray("alpha", Null.DIMENSION, new int[]{chiData.getLength()});
            }
            double[] chi = chiData.getData();
            AccumulatorVirialOverlapSingleAverage acc = dsvo.getAccumulators()[0];
            for (int i=0; i<chi.length; i++) {
                chi[i] = acc.getBennetBias(i);
            }
            return chiData;
        }

        public DataTag getTag() {
            return chiTag;
        }

        public IDataInfo getDataInfo() {
            if (chiInfo == null || chiInfo.getLength() != dsvo.getAccumulators()[0].getNBennetPoints()) {
                chiData = new DataFunction(new int[]{dsvo.getAccumulators()[0].getNBennetPoints()});
                alphaData = new DataDoubleArray(chiData.getLength());
                chiInfo = new DataInfoFunction("alpha", Null.DIMENSION, this);
            }
            return chiInfo;
        }

        public DataDoubleArray getIndependentData(int j) {
            double[] alpha = alphaData.getData();
            AccumulatorVirialOverlapSingleAverage acc = dsvo.getAccumulators()[0];
            for (int i=0; i<alpha.length; i++) {
                alpha[i] = acc.getBennetBias(i);
            }
            return alphaData;
        }

        public DataInfoDoubleArray getIndependentDataInfo(int i) {
            if (alphaInfo == null || alphaInfo.getLength() != dsvo.getAccumulators()[0].getNBennetPoints()) {
                alphaInfo = new DataInfoDoubleArray("alpha", Null.DIMENSION, new int[]{dsvo.getAccumulators()[0].getNBennetPoints()});
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
