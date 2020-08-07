/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.vle;

import java.util.HashMap;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.border.TitledBorder;

import etomica.action.IAction;
import etomica.action.activity.Controller;
import etomica.action.controller.Controller2;
import etomica.graphics.DeviceButton;
import etomica.graphics.DevicePlotPoints;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.units.Debye;
import etomica.units.Kelvin;
import etomica.math.function.Function;


public class B2Fit extends SimulationPanel {

	private final static String APP_NAME = "CO2 Potential Fitter";

	public DevicePlotPoints dPlot;

	public B2Fit() {
		super(APP_NAME);

        Controller c = new Controller(); // TODO
        c.controller2 = new Controller2();
        //
        // Function
        //
        FunctionB2LJQ functionB2LJQ = new FunctionB2LJQ();
        final FunctionMayerB2LJQ functionMayerB2LJQ = new FunctionMayerB2LJQ();

    	dPlot = new DevicePlotPoints(c, new String[] {"Q", "epsilon", "sigma"},
    			new Function[]{functionB2LJQ}, new String[] {"B2"}, false);
    	dPlot.setTableColumnNames(new String[]{"T (K)","B2 (mL/mol)"});
    	dPlot.setSliderTextboxesEditable(true);

    	functionB2LJQ.dPlot = dPlot;
//    	functionMayerB2LJQ.dPlot = dPlot;

        controlPanel.add(dPlot.controlGraphic());
		graphicsPanel.add(dPlot.graphic());
		footerPanel.add(dPlot.parameterGraphic());

        dPlot.setXScale(200, 500, 300);
        dPlot.getPlotSizeSlider(DevicePlotPoints.MIN_X).setValue(200);
        dPlot.getPlotSizeSlider(DevicePlotPoints.MAX_X).setValue(400);
        dPlot.setAutoScale(true);

        dPlot.setYScale(-100, 100, 0);
		
        dPlot.getSlider("Q").setPrecision(2);
        dPlot.getSlider("Q").setValue(5.0);
		dPlot.setParameterLimits("Q", 0, 10.0);
        dPlot.setParameterLimits("epsilon", 100.0, 300.0);
        dPlot.getSlider("epsilon").setValue(200);
        dPlot.getSlider("epsilon").setNMajor(2);
        dPlot.getSlider("sigma").setPrecision(2);
        dPlot.setParameterLimits("sigma", 2.0, 5.0);

        ((JPanel)dPlot.getSlider("Q").graphic().getParent()).setBorder(new TitledBorder(null, "Q (Debye-A)",
                TitledBorder.CENTER, TitledBorder.TOP));
        ((JPanel)dPlot.getSlider("sigma").graphic().getParent()).setBorder(new TitledBorder(null, "sigma (A)",
                TitledBorder.CENTER, TitledBorder.TOP));
        ((JPanel)dPlot.getSlider("epsilon").graphic().getParent()).setBorder(new TitledBorder(null, "epsilon (K)",
                TitledBorder.CENTER, TitledBorder.TOP));

        DeviceButton recalcButton = new DeviceButton(c);
        recalcButton.setAction(new IAction() {
            public void actionPerformed() {
                functionMayerB2LJQ.reset();
                dPlot.refresh();
            }
        });
        recalcButton.setLabel("Recalc Mayer");
//        controlPanel.add(recalcButton.graphic());

        dPlot.getDisplayPlot().getPlot().setXLabel("Temperature (K)");
        dPlot.getDisplayPlot().getPlot().setYLabel("B2 (mL/mol)");
    }

	/**
	 * Function using LJQB2 to evaluate B2 as a function of Q, epsilon, sigma
	 * and T.
	 */
	public static class FunctionB2LJQ implements Function {
	    public double f(double T) {
	        double Q = 0;
	        double epsilon = 1;
            double sigma = 1;
	        if (dPlot != null) {
	            Q = Debye.UNIT.toSim(dPlot.getParameterValue("Q"));
                epsilon = dPlot.getParameterValue("epsilon");
                sigma = dPlot.getParameterValue("sigma");
	        }
            T /= epsilon;
	        Q /= Math.sqrt(Kelvin.UNIT.toSim(epsilon))*Math.pow(sigma, 2.5);
	        return LJQB2.B2(T, Q)/1E24*etomica.util.Constants.AVOGADRO*Math.pow(sigma,3);
	    }
	    
	    public DevicePlotPoints dPlot;
	}

    /**
     * Function using Mayer sampling to evaluate B2 as a function of Q,
     * epsilon, sigma and T.
     */
    public static class FunctionMayerB2LJQ implements Function {
        public double f(double T) {
            // only evaluate B2 if user has hit the button
            if (!ready) return Double.NaN;
            double Q = 0;
            double epsilon = 1;
            double sigma = 1;
            if (dPlot != null) {
                Q = Debye.UNIT.toSim(dPlot.getParameterValue("Q"));
                epsilon = dPlot.getParameterValue("epsilon");
                sigma = dPlot.getParameterValue("sigma");
            }
            if (Double.isNaN(oldQ)) {
                oldQ = Q;
            }
            if (Double.isNaN(oldEpsilon)) {
                oldEpsilon = epsilon;
            }
            if (Double.isNaN(oldSigma)) {
                oldEpsilon = epsilon;
            }
            // try to use a cached value if Q, epsilon, sigma and T haven't changed
            if (Q == oldQ && epsilon == oldEpsilon && sigma == oldSigma) {
                Double B2 = B2Hash.get(T);
                if (B2 != null) {
                    return B2.doubleValue()*Math.pow(sigma,3);
                }
            }
            else {
                // user changed a parameter, but hasn't hit the recalc button
                return Double.NaN;
            }
            double unreducedT = T;
            T /= epsilon;
            Q /= Math.sqrt(Kelvin.UNIT.toSim(epsilon))*Math.pow(sigma, 2.5);
            double b2 = VirialLJQB2.calcB2(T, Q)/1E24*etomica.util.Constants.AVOGADRO;
            B2Hash.put(unreducedT, b2);
            return b2*Math.pow(sigma,3);
        }
        
        public void reset() {
            // forget cached B2, forget old parameter values
            B2Hash.clear();
            oldQ = Double.NaN;
            oldEpsilon = Double.NaN;
            oldSigma = Double.NaN;
            ready = true;
        }
        
        public DevicePlotPoints dPlot;
        protected HashMap<Double,Double> B2Hash = new HashMap<Double,Double>();
        protected double oldQ, oldEpsilon, oldSigma;
        public boolean ready = false;
    }

    public static void main(String[] args) {

        SimulationGraphic.initGraphics();
        B2Fit graph = new B2Fit();
        JFrame f = new JFrame();
        f.setSize(1000, 600);
        f.getContentPane().add(graph);
        f.setVisible(true);
        f.addWindowListener(SimulationGraphic.WINDOW_CLOSER);
        graph.dPlot.refresh();

    }

    public static class Applet extends javax.swing.JApplet {
	    public void init() {

            B2Fit graph = new B2Fit();
		    getContentPane().add(graph);
		    graph.dPlot.refresh();
	    }

        private static final long serialVersionUID = 1L;
    }

}
