package etomica.graphics;

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
	
	private static final double  DELAY_MULTIPLIER = 10;

    public DeviceDelaySlider(Controller cont, ActivityIntegrate ai) {

    	DelayModifier mod = new DelayModifier(ai);

    	delaySlider = new DeviceSlider(cont, mod);
    	delaySlider.setPrecision(1);
    	delaySlider.setMinimum(0);
    	delaySlider.setMaximum(10);
    	delaySlider.setNMajor(5);
    	delaySlider.setValue(0);

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
	    	return activityIntegrate.getSleepPeriod() / DELAY_MULTIPLIER;
	    }

	    public void setValue(double d) {
	    	activityIntegrate.setSleepPeriod((int)(d*DELAY_MULTIPLIER));
	    }

        public Dimension getDimension() {return Null.DIMENSION;}
        public String getLabel() {return "Sleep period";}

	}
}
