package etomica.multiharmonic;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;

import javax.swing.JPanel;

import etomica.action.Action;
import etomica.data.DataProcessorFunction;
import etomica.data.DataPump;
import etomica.data.DataSourceFunction;
import etomica.data.types.CastArrayToDoubleArray;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceTrioControllerButton;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxesCAE;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntervalActionAdapter;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.space1d.Vector1D;
import etomica.units.Dimension;
import etomica.util.Function;


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

        final DisplayPhase displayPhase = new DisplayPhase(sim.phase);
        sim.integrator.addListener(new IntervalActionAdapter(
                new Action() {
                    public void actionPerformed() {
                        displayPhase.repaint();
                    }
                    public String getLabel() {return "";}
                }));
        
        JPanel controlPanel = new JPanel(new GridBagLayout());
        controlPanel.add(control.graphic());
        
        panel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc2 = new GridBagConstraints();
        gbc2.gridx = 0;
        gbc2.gridy = java.awt.GridBagConstraints.RELATIVE;
        panel.add(controlPanel, gbc2);
        panel.add(displayPhase.graphic(), gbc2);
        
        DisplayBoxesCAE displayBox = new DisplayBoxesCAE();
        displayBox.setLabel("Free Energy");
        displayBox.setLabelType(DisplayBox.BORDER);
        displayBox.setAccumulator(sim.accumulator);
        panel.add(displayBox.graphic());
        
        DisplayPlot plot = new DisplayPlot();
        CastArrayToDoubleArray cast = new CastArrayToDoubleArray();
        DataProcessorFunction log = new DataProcessorFunction(Function.Log.INSTANCE);
        sim.history.setDataSink(cast);
        cast.setDataSink(log);
        log.setDataSink(plot.getDataTable());
        panel.add(plot.graphic());
        
        DeviceSlider x0Slider = new DeviceSlider(sim.controller);
        final DeviceSlider omegaSlider = new DeviceSlider(sim.controller);
        x0Slider.setPrecision(1);
        omegaSlider.setPrecision(1);
        Modifier x0Modifier = new Modifier() {
            public void setValue(double value) {
                sim.potentialB.setX0(new Vector1D(value));
            }
            public double getValue() {
                return sim.potentialB.getX0().x(0);
            }
            public String getLabel() {return "x0";}
            public Dimension getDimension() {return Dimension.LENGTH;}
        };
        x0Slider.setModifier(x0Modifier);
        x0Slider.setMinimum(0.0);
        x0Slider.setMaximum(5.0);
        omegaSlider.setModifier(new ModifierGeneral(sim.potentialB, "springConstant"));
        omegaSlider.setMinimum(1.0);
        omegaSlider.setMaximum(100.0);
        JPanel sliderPanel = new JPanel(new GridLayout(0,2));
        sliderPanel.add(x0Slider.graphic());
        sliderPanel.add(omegaSlider.graphic());
        controlPanel.add(sliderPanel, gbc2);
        omegaSlider.getSlider().setBorder(new javax.swing.border.TitledBorder("omegaB"));
        x0Slider.getSlider().setBorder(new javax.swing.border.TitledBorder("x0"));
        
        Function deltaF = new Function() {
            public double f(double x) {
                System.out.println(omegaSlider.getValue());
                return sim.phase.atomCount() * Math.log(omegaSlider.getValue());
            }
            public double dfdx(double x) {
                return 0.0;
            }
            public double inverse(double f) {
                return 0.0;
            }
        };
        DataSourceFunction exact = new DataSourceFunction(deltaF);
        DataPump exactPump = new DataPump(exact, plot.getDataTable());
        omegaSlider.setPostAction(exactPump);
        plot.getDataTable().setUpdatingOnAnyChange(true);


    }
    
    public static void main(String[] args) {
        MultiharmonicGraphic simGraphic = new MultiharmonicGraphic();
        SimulationGraphic.makeAndDisplayFrame(simGraphic.panel);
    }


    JPanel panel;

}
