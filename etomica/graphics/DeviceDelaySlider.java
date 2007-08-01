package etomica.graphics;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.TitledBorder;

import etomica.action.activity.Controller;
import etomica.action.activity.ActivityIntegrate;
import etomica.modifier.Modifier;
import etomica.units.Dimension;
import etomica.units.Null;

public class DeviceDelaySlider {

	private DeviceSlider      delaySlider;
	private JPanel            delayPanel;
	
    //DELAY_EXPONENT affects how sharply the delay increases as slider is moved from zero -- 
    //a larger value pushes increase off to larger slider values; 1.0 is a linear increase
    //DELAY_MULTIPLIER is set such that sleep period is 100 when slider is at its maximum value of 10
    private static final double  DELAY_EXPONENT = 2.0;
    private static final double  DELAY_MULTIPLIER = 100.0 / Math.pow(10.0,DELAY_EXPONENT);

    public DeviceDelaySlider(Controller cont, ActivityIntegrate ai) {

    	DelayModifier mod = new DelayModifier(ai);

    	delaySlider = new DeviceSlider(cont, mod);
    	delaySlider.setShowValues(false);
    	delaySlider.setPrecision(1);
    	delaySlider.setMinimum(0);
    	delaySlider.setMaximum(10);
    	delaySlider.setValue(0);
    	delaySlider.setNMajor(0);
	    java.util.Hashtable scaleLabels = new java.util.Hashtable();
	    scaleLabels.put(new Integer(0), new JLabel( "fast", JLabel.CENTER ));
	    // slow is 100 : need to know details of DeviceSlider to understand why.
	    scaleLabels.put(new Integer(100), new JLabel( "slow", JLabel.CENTER ));
	    delaySlider.getSlider().setLabelTable(scaleLabels);

    	delayPanel = new JPanel();
    	delayPanel.setBorder(new TitledBorder(null, "Simulation Delay", TitledBorder.CENTER, TitledBorder.TOP));

    	delayPanel.add(delaySlider.graphic());

    }

	public JPanel graphic() {
		return delayPanel;
	}

	private class DelayModifier implements Modifier {

        ActivityIntegrate activityIntegrate;

        public DelayModifier(ActivityIntegrate ai) {
    	    activityIntegrate = ai;
        }

	    public double getValue() {
	    	return Math.pow(activityIntegrate.getSleepPeriod() / DELAY_MULTIPLIER, 1.0/DELAY_EXPONENT);
	    }

	    public void setValue(double d) {
	    	activityIntegrate.setSleepPeriod((int)(Math.pow(d, DELAY_EXPONENT)*DELAY_MULTIPLIER));
	    }

        public Dimension getDimension() {return Null.DIMENSION;}
        public String getLabel() {return "Sleep period";}

	}
}
