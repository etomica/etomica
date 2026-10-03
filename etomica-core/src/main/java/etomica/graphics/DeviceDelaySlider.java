/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;


import etomica.action.controller.Controller;
import etomica.modifier.Modifier;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Null;

import javax.swing.*;
import javax.swing.border.TitledBorder;

public class DeviceDelaySlider {

	private final JPanel delayPanel;
    private final Controller controller;

    //DELAY_EXPONENT affects how sharply the delay increases as slider is moved from zero --
    //a larger value pushes increase off to larger slider values; 1.0 is a linear increase
    //DELAY_MULTIPLIER is set such that sleep period is 100 when slider is at its maximum value of 10
    protected double  delayExponent;
    protected int maxSleep;

	//TODO develop means to optionally reduce integrator time step when delay is added, for smoother or more complete animation (e.g., advance spheres closer to collision point)

    public DeviceDelaySlider(Controller cont) {
		this.controller = cont;
        delayExponent = 2.0;
        maxSleep = 100;
    	DelayModifier mod = new DelayModifier();

		DeviceSlider delaySlider = new DeviceSlider(cont, mod);
    	delaySlider.setShowValues(false);
    	delaySlider.setPrecision(1);
    	delaySlider.setMinimum(0);
    	delaySlider.setMaximum(10);
    	delaySlider.setValue(0);
    	delaySlider.setNMajor(0);
	    java.util.Hashtable<Integer,JLabel> scaleLabels = new java.util.Hashtable<Integer,JLabel>();
	    scaleLabels.put(0, new JLabel( "fast", JLabel.CENTER ));
	    // slow is 100 : need to know details of DeviceSlider to understand why.
	    scaleLabels.put(100, new JLabel( "slow", JLabel.CENTER ));
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

	public JPanel graphic() {
		return delayPanel;
	}

	private class DelayModifier implements Modifier {

	    public double getValue() {
	        double delayMultiplier = maxSleep / Math.pow(10.0,delayExponent);
	    	return Math.pow(controller.getSleepPeriod() / delayMultiplier, 1.0/delayExponent);
	    }

	    public void setValue(double d) {
			double delayMultiplier = maxSleep / Math.pow(10.0, delayExponent);
			double sleep = Math.pow(d, delayExponent) * delayMultiplier;
			controller.setSleepPeriod(sleep);
		}

        public Dimension getDimension() {return Null.DIMENSION;}
        public String getLabel() {return "Sleep period";}

	}
}
