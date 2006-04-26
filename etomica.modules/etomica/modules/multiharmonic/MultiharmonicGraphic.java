package etomica.modules.multiharmonic;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;

import javax.swing.JPanel;

import etomica.action.Action;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorHistory;
import etomica.data.DataProcessorFunction;
import etomica.data.DataPump;
import etomica.data.DataSourceFunction;
import etomica.data.DataSourceScalar;
import etomica.data.types.CastGroupToDoubleArray;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceTrioControllerButton;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntervalActionAdapter;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.space1d.Vector1D;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Length;
import etomica.units.Null;
import etomica.units.Pixel;
import etomica.util.Function;
import etomica.util.HistoryCollapsing;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Oct 28, 2005 by kofke
 */
public class MultiharmonicGraphic {

    /**
     * 
     */
    public MultiharmonicGraphic() {
        super();
        final Multiharmonic sim = new Multiharmonic();

        sim.register(sim.integrator);
        DeviceTrioControllerButton control = new DeviceTrioControllerButton(sim);

        final DisplayPhase displayPhase = new DisplayPhase(sim.phase,sim.getDefaults().pixelUnit);
        sim.integrator.addListener(new IntervalActionAdapter(
                new Action() {
                    public void actionPerformed() {
                        displayPhase.repaint();
                    }
                    public String getLabel() {return "";}
                }));
        displayPhase.setPixelUnit(new Pixel(400/sim.phase.getBoundary().getDimensions().x(0)));

        
 //       panel.add(displayPhase.graphic(), gbc2);
        
//        DisplayBoxesCAE displayBox = new DisplayBoxesCAE();
//        displayBox.setLabel("Free Energy");
//        displayBox.setLabelType(DisplayBox.BORDER);
//        displayBox.setAccumulator(sim.accumulator);
//        panel.add(displayBox.graphic());
        
        DisplayPlot plot = new DisplayPlot();
        DataProcessorFunction log = new DataProcessorFunction(new Function() {
          public double f(double x) {
              return -Math.log(x);
          }
          public double dfdx(double x) {
              return 0.0;
          }
          public double inverse(double f) {
              return 0.0;
          }});
        sim.accumulator.addDataSink(log, new AccumulatorAverage.StatType[] {AccumulatorAverage.StatType.AVERAGE});
        AccumulatorHistory history = new AccumulatorHistory(HistoryCollapsing.FACTORY);
        history.setTimeDataSource(sim.timeCounter);
        log.setDataSink(history);
        history.setDataSink(plot.getDataSet());
        
        DisplayPlot energyPlot = new DisplayPlot();
        sim.historyEnergy.setTimeDataSource(sim.timeCounter);
        sim.historyEnergy.setDataSink(energyPlot.getDataSet());
        
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
                return sim.potentialB.getX0().x(0);
            }
            public String getLabel() {return "x0";}
            public Dimension getDimension() {return Length.DIMENSION;}
        };
        x0Slider.setModifier(x0Modifier);
        x0Slider.setMinimum(1.0);
        x0Slider.setMaximum(3.0);
        x0Slider.setValue(1.0);
        omegaASlider.setModifier(new ModifierGeneral(sim.potentialA, "springConstant"));
        omegaASlider.setMinimum(0.1);
        omegaASlider.setMaximum(50.0);
        omegaASlider.setValue(1.0);
        omegaBSlider.setModifier(new ModifierGeneral(sim.potentialB, "springConstant"));
        omegaBSlider.setMinimum(0.1);
        omegaBSlider.setMaximum(10.0);
        omegaBSlider.setValue(1.0);

        omegaASlider.getSlider().setBorder(new javax.swing.border.TitledBorder("omegaA"));
        omegaBSlider.getSlider().setBorder(new javax.swing.border.TitledBorder("omegaB"));
        x0Slider.getSlider().setBorder(new javax.swing.border.TitledBorder("x0"));
        
        DataSourceScalar delta = new DataSourceScalar("exact",Energy.DIMENSION) {
            public double getDataAsScalar() {
                return 0.5*sim.phase.atomCount() * Math.log(omegaBSlider.getValue()/omegaASlider.getValue());
            }
        };
        DataSourceScalar uAvg = new DataSourceScalar("exact",Energy.DIMENSION) {
            public double getDataAsScalar() {
                return sim.phase.atomCount();
            }
        };
        
        AccumulatorHistory deltaHistory = new AccumulatorHistory(HistoryCollapsing.FACTORY,sim.historyEnergy.getDataLength());
        DataPump exactPump = new DataPump(delta, deltaHistory);
        deltaHistory.setDataSink(plot.getDataSet());
        IntervalActionAdapter adapter = new IntervalActionAdapter(exactPump);
        adapter.setActionInterval(sim.accumulator.getBlockSize());
        sim.integrator.addListener(adapter);
        sim.register(delta, exactPump);
        deltaHistory.setTimeDataSource(sim.timeCounter);
        
        AccumulatorHistory uAvgHistory = new AccumulatorHistory(HistoryCollapsing.FACTORY,sim.historyEnergy.getDataLength());
        DataPump uPump = new DataPump(uAvg, uAvgHistory);
        uAvgHistory.setDataSink(energyPlot.getDataSet());
        adapter = new IntervalActionAdapter(uPump);
        adapter.setActionInterval(sim.accumulatorEnergy.getBlockSize());
        sim.integrator.addListener(adapter);
        sim.register(uAvg, uPump);
        uAvgHistory.setTimeDataSource(sim.timeCounter);
        
        plot.getDataSet().setUpdatingOnAnyChange(true);
        energyPlot.getDataSet().setUpdatingOnAnyChange(true);
        plot.getPlot().setTitle("Free energy difference");
        energyPlot.getPlot().setTitle("Average energy");

        final DisplayPlot uPlot = new DisplayPlot();
        final double yMax = 2.0;
        uPlot.getPlot().setYRange(0.0, yMax);
        
        Function fUA = new Function() {
            public double f(double x) {
                double x0 = sim.potentialA.getX0().x(0);
                return 0.5*sim.potentialA.getSpringConstant()*(x - x0)*(x - x0);
            }
            public double dfdx(double x) {
                return 0.0;
            }
            public double inverse(double f) {
                return 0.0;
            }
        };
        Function fUB = new Function() {
            public double f(double x) {
                double x0 = sim.potentialB.getX0().x(0);
                return 0.5*sim.potentialB.getSpringConstant()*(x - x0)*(x - x0);
            }
            public double dfdx(double x) {
                return 0.0;
            }
            public double inverse(double f) {
                return 0.0;
            }
        };

        final DataSourceFunction uA = new DataSourceFunction("A",Null.DIMENSION,fUA,100,"x",Length.DIMENSION);
        final DataSourceFunction uB = new DataSourceFunction("B",Null.DIMENSION,fUB,100,"x",Length.DIMENSION);
        uA.getXSource().setXMax(sim.phase.getBoundary().getDimensions().x(0));
        uB.getXSource().setXMax(sim.phase.getBoundary().getDimensions().x(0));
        uAPump = new DataPump(uA, uPlot.getDataSet());
        uBPump = new DataPump(uB, uPlot.getDataSet());
        Action uUpdate = new Action() {
            public void actionPerformed() {
                uA.update();
                uB.update();
                uAPump.actionPerformed();
                uBPump.actionPerformed();
            }
            public String getLabel() {return "";}
        };
        omegaASlider.setPostAction(uUpdate);
        omegaBSlider.setPostAction(uUpdate);
        x0Slider.setPostAction(uUpdate);
        
        uPlot.getDataSet().setUpdatingOnAnyChange(true);

        //Lay out components
        //main panel
        panel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc2 = new GridBagConstraints();
        gbc2.gridx = 0;
        gbc2.gridy = 0;
        
        //controls -- start/pause and sliders
        JPanel controlPanel = new JPanel(new GridBagLayout());
        panel.add(controlPanel, gbc2);
        controlPanel.add(control.graphic(), gbc2);
        JPanel sliderPanel = new JPanel(new GridLayout(3,1));
        sliderPanel.add(x0Slider.graphic());
        sliderPanel.add(omegaASlider.graphic());
        sliderPanel.add(omegaBSlider.graphic());
        gbc2.gridx = 0;
        gbc2.gridy = 1;
        controlPanel.add(sliderPanel, gbc2);
        
        //energy plot
        energyPlot.setSize(300,200);
        panel.add(energyPlot.graphic(),gbc2);
        
        //plot of potential and display of phase
        JPanel phasePanel = new JPanel(new GridBagLayout());
        gbc2.gridx = 0;
        gbc2.gridy = 0;
        phasePanel.add(uPlot.graphic(), gbc2);
        gbc2.gridx = 0;
        gbc2.gridy = 1;
        phasePanel.add(displayPhase.graphic(), gbc2);
        gbc2.gridx = 1;
        gbc2.gridy = 0;
        panel.add(phasePanel, gbc2);
        
        //plot of running average
        gbc2.gridx = 1;
        gbc2.gridy = 1;
        panel.add(plot.graphic(), gbc2);

        uUpdate.actionPerformed();
        
    }
    
    /**
     * This pulls the data from the DataSourceFunctions and pushes them to the
     * plot.  Doing this before was ineffective because the plot was not visible
     * and DispalyPlot refuses to update in that situation. 
     */
    void initUPlot() {
        uAPump.actionPerformed();
        uBPump.actionPerformed();
    }
    
    public static void main(String[] args) {
        MultiharmonicGraphic simGraphic = new MultiharmonicGraphic();
        SimulationGraphic.makeAndDisplayFrame(simGraphic.panel);
        simGraphic.initUPlot();
    }

    public static class Applet extends javax.swing.JApplet {

        public void init() {
            MultiharmonicGraphic simGraphic = new MultiharmonicGraphic();
            getContentPane().add(simGraphic.panel);
            simGraphic.initUPlot();
        }
    }

    JPanel panel;
    DataPump uAPump;
    DataPump uBPump;
}
