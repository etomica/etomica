/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.TitledBorder;

import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.IController;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorMD;
import etomica.modifier.Modifier;
import etomica.units.Dimension;
import etomica.units.Null;

public class DeviceDelaySlider {

	private DeviceSlider      delaySlider;
	private JPanel            delayPanel;
	protected double threshold = Double.POSITIVE_INFINITY;
	protected double nominalTimeStep;
	protected final ActivityIntegrate ai;
	
    //DELAY_EXPONENT affects how sharply the delay increases as slider is moved from zero -- 
    //a larger value pushes increase off to larger slider values; 1.0 is a linear increase
    //DELAY_MULTIPLIER is set such that sleep period is 100 when slider is at its maximum value of 10
    protected double  delayExponent;
    protected int maxSleep;

    public DeviceDelaySlider(IController cont, ActivityIntegrate ai) {
        delayExponent = 2.0;
        maxSleep = 100;
        this.ai = ai;
    	DelayModifier mod = new DelayModifier();

    	delaySlider = new DeviceSlider(cont, mod);
    	delaySlider.setShowValues(false);
    	delaySlider.setPrecision(1);
    	delaySlider.setMinimum(0);
    	delaySlider.setMaximum(10);
    	delaySlider.setValue(0);
    	delaySlider.setNMajor(0);
	    java.util.Hashtable<Integer,JLabel> scaleLabels = new java.util.Hashtable<Integer,JLabel>();
	    scaleLabels.put(new Integer(0), new JLabel( "fast", JLabel.CENTER ));
	    // slow is 100 : need to know details of DeviceSlider to understand why.
	    scaleLabels.put(new Integer(100), new JLabel( "slow", JLabel.CENTER ));
	    delaySlider.getSlider().setLabelTable(scaleLabels);

    	delayPanel = new JPanel();
    	delayPanel.setBorder(new TitledBorder(null, "Simulation Delay", TitledBorder.CENTER, TitledBorder.TOP));

    	delayPanel.add(delaySlider.graphic());
    }
    
    public void setDelayExponent(double newDelayExponent) {
        delayExponent = newDelayExponent;
    }
    
    public void setMaxSleep(int newMaxSleep) {
        maxSleep = newMaxSleep;
    }

    public double getDelayExponent() {
        return delayExponent;
    }
    
    public int getMaxSleep() {
        return maxSleep;
    }
    
    /**
     * You play with fire!  This allows the slider to change not increase
     * the activity integrate sleep period, but also to decrease the integrator
     * timestep.  This will almost certainly mess up data collection (averages,
     * etc. need to be reset when moving the slider whenever the slider changes
     * the timestep).
     * 
     * @param t The maximum sleep period that will actually be invoked.  If the
     * desired sleep period is greater than t, then the sleep period will be
     * equal to t and the integrator timestep will be lowered so that the amount
     * slept per simulation time is as desired.
     */
    public void setThreshold(int t) {
        Integrator integrator = ai.getIntegrator();
        if (integrator instanceof IntegratorMD) {
            nominalTimeStep = ((IntegratorMD)integrator).getTimeStep();
        }
        else {
            throw new RuntimeException("integrator is not MD, cannot set threshold");
        }
        if (t < 0.1) {
            throw new RuntimeException("threshold must be at least 0.1");
        }
        threshold = t;
    }
    
	public JPanel graphic() {
		return delayPanel;
	}

	private class DelayModifier implements Modifier {

	    public double getValue() {
	        double delayMultiplier = maxSleep / Math.pow(10.0,delayExponent);
	    	return Math.pow(ai.getSleepPeriod() / delayMultiplier, 1.0/delayExponent);
	    }

	    public void setValue(double d) {
	        double delayMultiplier = maxSleep / Math.pow(10.0,delayExponent);
            double sleep = Math.pow(d, delayExponent)*delayMultiplier;
            double tStep = nominalTimeStep;
            if (sleep > threshold) {
                tStep /= (sleep/threshold);
                sleep = threshold;
            }
            ai.setSleepPeriod((int)sleep);
            if (nominalTimeStep>0) ((IntegratorMD)ai.getIntegrator()).setTimeStep(tStep);
	    }

        public Dimension getDimension() {return Null.DIMENSION;}
        public String getLabel() {return "Sleep period";}

	}
}
