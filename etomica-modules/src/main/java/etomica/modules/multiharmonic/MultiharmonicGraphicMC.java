/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic;

import etomica.action.IAction;
import etomica.atom.DiameterHashByType;
import etomica.data.*;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.histogram.HistogramCollapsing;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.graphics.*;
import etomica.listener.IntegratorListenerAction;
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
        super(simulation, GRAPHIC_ONLY, APP_NAME, REPAINT_INTERVAL, _space, simulation.getController());
        this.sim = simulation;
        final DisplayBox displayBox = getDisplayBox(sim.box);
        remove(displayBox);
        final Device controllerButtons = getController();
        getPanel().graphicsPanel.remove(controllerButtons.graphic());
        getPanel().footerPanel.add(controllerButtons.graphic());
        getPanel().graphicsPanel.setLayout(new GridBagLayout());
        
        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();
        dataStreamPumps.add(simulation.dataPump);
        dataStreamPumps.add(simulation.dataPumpEnergy);

        displayBox.setPixelUnit(new Pixel(350/sim.box.getBoundary().getBoxSize().getX(0)));
        ((DiameterHashByType)displayBox.getDiameterHash()).setDiameter(sim.species.getLeafType(), 0.02);

        final DisplayPlot plot = new DisplayPlot();
        DataProcessorFunction log = new DataProcessorFunction(new Function() {
            public double f(double x) {return -Math.log(x);}
        });
        sim.accumulator.addDataSink(log, new AccumulatorAverage.StatType[] {sim.accumulator.AVERAGE});
        sim.accumulator.setPushInterval(1);
        AccumulatorHistory history = new AccumulatorHistory(new HistoryCollapsingDiscard(102, 3));
        history.setTimeDataSource(sim.stepCounter);
        log.setDataSink(history);
        history.setPushInterval(1000);
        history.setDataSink(plot.getDataSet().makeDataSink());
        plot.setLegend(new DataTag[]{history.getTag()}, "measured");
        
        final DisplayPlot energyPlot = new DisplayPlot();
        sim.historyEnergy.setTimeDataSource(sim.stepCounter);
        sim.historyEnergy.setDataSink(energyPlot.getDataSet().makeDataSink());
        energyPlot.setLegend(new DataTag[]{sim.historyEnergy.getTag()}, "measured");
        sim.historyEnergy.setPushInterval(10);
        
        DataProcessorFunction deltaU = new DataProcessorFunction(Function.Log.INSTANCE);
        sim.accumulator.addDataSink(deltaU, new StatType[]{sim.accumulator.MOST_RECENT});
        AccumulatorHistogram duHistogram = new AccumulatorHistogram(new HistogramCollapsing());
        deltaU.setDataSink(duHistogram);
        duHistogram.setPushInterval(1000);
        final DisplayPlot duPlot = new DisplayPlot();
        duHistogram.setDataSink(duPlot.getDataSet().makeDataSink());
        duPlot.setDoLegend(false);
        
        DeviceSlider x0Slider = new DeviceSlider(sim.controller);
        final DeviceSlider omegaASlider = new DeviceSlider(sim.controller);
        final DeviceSlider omegaBSlider = new DeviceSlider(sim.controller);
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
                return 0.5*sim.box.getLeafList().getAtomCount() * Math.log(omegaBSlider.getValue()/omegaASlider.getValue());
            }
        };
        DataSourceScalar uAvg = new DataSourceScalar("exact",Energy.DIMENSION) {
            public double getDataAsScalar() {
                return 0.5*sim.box.getLeafList().getAtomCount();
            }
        };
        
        int historyLength = sim.historyEnergy.getHistory().getHistoryLength();
        int nCollapseBins = ((HistoryCollapsingDiscard)sim.historyEnergy.getHistory()).getNumCollapseBins();
        AccumulatorHistory deltaHistory = new AccumulatorHistory(new HistoryCollapsingDiscard(historyLength, nCollapseBins));
        deltaHistory.setPushInterval(10000);
        DataPump exactPump = new DataPump(delta, deltaHistory);
        deltaHistory.setDataSink(plot.getDataSet().makeDataSink());
        IntegratorListenerAction exactPumpListener = new IntegratorListenerAction(exactPump);
        sim.integrator.getEventManager().addListener(exactPumpListener);
        exactPumpListener.setInterval((int)sim.accumulator.getBlockSize());
        dataStreamPumps.add(exactPump);
        deltaHistory.setTimeDataSource(sim.stepCounter);
        
        AccumulatorHistory uAvgHistory = new AccumulatorHistory(new HistoryCollapsingDiscard(historyLength, nCollapseBins));
        uAvgHistory.setPushInterval(10000);
        DataPump uPump = new DataPump(uAvg, uAvgHistory);
        uAvgHistory.setDataSink(energyPlot.getDataSet().makeDataSink());
        IntegratorListenerAction uPumpListener = new IntegratorListenerAction(uPump);
        sim.integrator.getEventManager().addListener(uPumpListener);
        uPumpListener.setInterval((int)sim.accumulatorEnergy.getBlockSize());
        dataStreamPumps.add(uPump);
        uAvgHistory.setTimeDataSource(sim.stepCounter);
        
        plot.getDataSet().setUpdatingOnAnyChange(true);
        energyPlot.getDataSet().setUpdatingOnAnyChange(true);
        plot.getPlot().setTitle("Free energy difference");
        energyPlot.getPlot().setTitle("Average energy");
        
        duPlot.getPlot().setTitle("Energy Difference (A-B)");

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
        uA.getXSource().setXMin(-sim.box.getBoundary().getBoxSize().getX(0));
        uB.getXSource().setXMin(-sim.box.getBoundary().getBoxSize().getX(0));
        uA.getXSource().setXMax(sim.box.getBoundary().getBoxSize().getX(0));
        uB.getXSource().setXMax(sim.box.getBoundary().getBoxSize().getX(0));
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

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();
        GridBagConstraints horizGBC = SimulationPanel.getHorizGBC();

        //controls -- start/pause and sliders
        JTabbedPane sliderPanel = new JTabbedPane();
        sliderPanel.add(x0Slider.graphic(), "x0");
        sliderPanel.add(omegaASlider.graphic(), "omegaA");
        sliderPanel.add(omegaBSlider.graphic(), "omegaB");
        
        energyPlot.setSize(350,250);
        duPlot.setSize(350,250);
        
        JPanel topPanel = new JPanel(new GridBagLayout());
        JPanel displayPanel = new JPanel(new GridBagLayout());
        topPanel.add(sliderPanel);
        topPanel.add(nSlider.graphic(), vertGBC);
        displayPanel.add(uPlot.graphic(), vertGBC);
        displayPanel.add(displayBox.graphic(), vertGBC);
        GridBagConstraints mygbc = new GridBagConstraints();
        mygbc.gridx = 1;
        mygbc.gridy = 0;
        mygbc.gridheight = 2;
        topPanel.add(displayPanel, mygbc);
        
        plot.setSize(350, 250);
        uPlot.setSize(450, 250);

        getPanel().graphicsPanel.add(topPanel);

        JTabbedPane plotTabs = new JTabbedPane();
        JPanel tab1 = new JPanel(new GridBagLayout());
        plotTabs.add(tab1, "Energy");
        tab1.add(energyPlot.graphic());
        tab1.add(duPlot.graphic(), horizGBC);

        getPanel().graphicsPanel.add(plotTabs, vertGBC);

        getController().getReinitButton().setPostAction(new IAction() {
        	public void actionPerformed() {
                displayBox.repaint();
                energyPlot.getPlot().repaint();
                plot.getPlot().repaint();
        	}
        });

        uUpdate.actionPerformed();

        AccumulatorAverageCollapsingLog accumulatorEnergy = new AccumulatorAverageCollapsingLog();
        DataPump dataPumpEnergy = new DataPump(sim.meter, accumulatorEnergy);
        dataStreamPumps.add(dataPumpEnergy);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(dataPumpEnergy));

        DataFork feFork = new DataFork();
        DataProcessorFunction negative = new DataProcessorFunction(new Function() {
            public double f(double x) {
                return -x;
            }
        });
        accumulatorEnergy.setDataSink(negative);
        negative.setDataSink(feFork);
        DisplayPlot funPlot = new DisplayPlot();
        funPlot.getPlot().setTitle("Free Energy Difference Convergence");
        funPlot.setDoLegend(false);
        funPlot.setSize(350, 250);
        feFork.addDataSink(funPlot.getDataSet().makeDataSink());
        funPlot.getPlot().setXLog(true);
        
        DataProcessorBounds dpUpper = new DataProcessorBounds();
        dpUpper.setAccumulator(accumulatorEnergy);
        dpUpper.setIsUpperBound(true);
        feFork.addDataSink(dpUpper);
        dpUpper.setDataSink(funPlot.getDataSet().makeDataSink());
        funPlot.setLegend(new DataTag[]{dpUpper.getTag()}, "FE upper bound");
        DataProcessorBounds dpLower = new DataProcessorBounds();
        dpLower.setAccumulator(accumulatorEnergy);
        dpLower.setIsUpperBound(false);
        feFork.addDataSink(dpLower);
        dpLower.setDataSink(funPlot.getDataSet().makeDataSink());
        funPlot.setLegend(new DataTag[]{dpLower.getTag()}, "FE lower bound");
        
        JPanel tab2 = new JPanel(new GridBagLayout());
        plotTabs.add(tab2, "Free Energy");
        tab2.add(funPlot.graphic());
        tab2.add(plot.graphic(), horizGBC);
        
        DataProcessorVar dpstdev = new DataProcessorVar();
        dpstdev.setAccumulator(accumulatorEnergy);
        feFork.addDataSink(dpstdev);
        DataFork stdevFork = new DataFork();
        dpstdev.setDataSink(stdevFork);
        DataProcessorFunction lnstdev = new DataProcessorFunction(new Function.Log());
        stdevFork.addDataSink(lnstdev);

        DataProcessorDyDLnx dAdlnx = new DataProcessorDyDLnx();
        feFork.addDataSink(dAdlnx);
        
        DataProcessorFunction ndAdlnx = new DataProcessorFunction(new Function() {
            public double f(double x) {return -x;}
        });
        dAdlnx.setDataSink(ndAdlnx);
        DataFork ndAdlnxFork = new DataFork();
        ndAdlnx.setDataSink(ndAdlnxFork);
        DataProcessorFunction lndAdlnx = new DataProcessorFunction(new Function.Log());
        ndAdlnxFork.addDataSink(lndAdlnx);
        DataProcessorDyDLnx dlndAdlnx = new DataProcessorDyDLnx();
        lndAdlnx.setDataSink(dlndAdlnx);

        DisplayPlot dfunPlot = new DisplayPlot();
        dfunPlot.getPlot().setTitle("d(ln(x))/dlnN");
        dfunPlot.setSize(350, 250);
        dlndAdlnx.setDataSink(dfunPlot.getDataSet().makeDataSink());
        dfunPlot.setLegend(new DataTag[]{dlndAdlnx.getTag()}, "x=-dA/dlnN");
        dfunPlot.getPlot().setXLog(true);
        dfunPlot.setDoLegend(true);

        DataProcessorBias dpBias = new DataProcessorBias(delta);
        feFork.addDataSink(dpBias);
        DataFork biasFork = new DataFork();
        dpBias.setDataSink(biasFork);
        DataProcessorFunction lnb = new DataProcessorFunction(new Function.Log());
        biasFork.addDataSink(lnb);
        
        DataProcessorDyDLnx dlnbDlnx = new DataProcessorDyDLnx();
        lnb.setDataSink(dlnbDlnx);
        
        dlnbDlnx.setDataSink(dfunPlot.getDataSet().makeDataSink());
        dfunPlot.setLegend(new DataTag[]{dlnbDlnx.getTag()}, "x=b");
        
        DataProcessorDyDLnx dlnstdevDlnx = new DataProcessorDyDLnx();
        lnstdev.setDataSink(dlnstdevDlnx);
        dlnstdevDlnx.setDataSink(dfunPlot.getDataSet().makeDataSink());
        dfunPlot.setLegend(new DataTag[]{dlnstdevDlnx.getTag()}, "x=var_A");

        DisplayPlot lfunPlot = new DisplayPlot();
        lfunPlot.setSize(300, 250);
        lfunPlot.getPlot().setTitle("x");
        ndAdlnxFork.addDataSink(lfunPlot.getDataSet().makeDataSink());
        lfunPlot.setLegend(new DataTag[]{ndAdlnxFork.getTag()}, "x=-dA/dlnN");
        lfunPlot.getPlot().setXLog(true);
        lfunPlot.getPlot().setYLog(true);

        stdevFork.addDataSink(lfunPlot.getDataSet().makeDataSink());
        lfunPlot.setLegend(new DataTag[]{stdevFork.getTag()}, "x=var_A");
        
        biasFork.addDataSink(lfunPlot.getDataSet().makeDataSink());
        lfunPlot.setLegend(new DataTag[]{biasFork.getTag()}, "x=b");

        JPanel tab3 = new JPanel(new GridBagLayout());
        plotTabs.add(tab3, "Bias");
        tab3.add(lfunPlot.graphic());
        tab3.add(dfunPlot.graphic(), horizGBC);

        sim.activityIntegrate.setSleepPeriod(0);
    }

    public static void main(String[] args) {
        final MultiharmonicMC sim = new MultiharmonicMC();
        MultiharmonicGraphicMC simGraphic = new MultiharmonicGraphicMC(sim, sim.getSpace());
        SimulationGraphic.makeAndDisplayFrame(simGraphic.getPanel(), APP_NAME);
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
