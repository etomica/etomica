/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import java.awt.Component;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;

import javax.swing.JPanel;
import javax.swing.border.TitledBorder;

import etomica.action.IAction;
import etomica.action.activity.Controller;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.units.Unit;


public class DeviceThermoSlider extends Device {

	protected final JPanel        temperaturePanel;  // main panel for thermo device PRIVATE
	protected final DeviceSlider  temperatureSlider; // Do not make make accessible
    protected final DeviceButtonGroup thermalButtons;

	protected static final int DEFAULT_MIN_TEMPERATURE = 0;
	protected static final int DEFAULT_MAX_TEMPERATURE = 300;

	public DeviceThermoSlider(Controller cont, final IntegratorBox integrator) {
        //adiabatic/isothermal radio button
	    thermalButtons = new DeviceButtonGroup(cont, 2);
	    thermalButtons.addButton("Adiabatic", new IAction() {
	        public void actionPerformed() {integrator.setIsothermal(false);}
	    });
        thermalButtons.addButton("Isothermal", new IAction() {
            public void actionPerformed() {integrator.setIsothermal(true);}
        });
        thermalButtons.setPostAction(new IAction() {
            public void actionPerformed() {
                configureSliderAccessibility();
            }
        });

        //temperature selector
        temperatureSlider = new DeviceSlider(controller);
        temperatureSlider.setShowValues(true);
        temperatureSlider.setEditValues(true);
        temperatureSlider.setMinimum(DEFAULT_MIN_TEMPERATURE);
        temperatureSlider.setMaximum(DEFAULT_MAX_TEMPERATURE);
        temperatureSlider.setNMajor(4);
        temperatureSlider.setValue(300);
        temperatureSlider.getSlider().setEnabled(false);
        temperatureSlider.getTextField().setEnabled(false);

        setController(cont);

        temperaturePanel = new JPanel(new GridBagLayout());
        temperaturePanel.setBorder(new TitledBorder(null, "Set Temperature", TitledBorder.CENTER, TitledBorder.TOP));
        GridBagConstraints gbc1 = new GridBagConstraints();
        gbc1.gridx = 0;  gbc1.gridy = 0;
        temperaturePanel.add(thermalButtons.graphic(), gbc1);
        gbc1.gridx = 0;  gbc1.gridy = 1;
        temperaturePanel.add(temperatureSlider.graphic(),gbc1);
        
        temperatureSlider.setModifier(new ModifierGeneral(integrator, "temperature"));
        if (integrator.isIsothermal()) {
            setIsothermal();
        }
        else {
            setAdiabatic();
        }
    }

	public void setIsothermalButtonsVisibility(boolean doShowIsothermalButtons) {
	    thermalButtons.graphic().setVisible(doShowIsothermalButtons);
	}

	public boolean getIsothermalButtonsVisibility() {
	    return thermalButtons.graphic().isVisible();
	}

	/**
	 * Set the Isothermal button to its selected state.
	 */
	public void setIsothermal() {
	    thermalButtons.setSelected("Isothermal");
	}

	/**
	 * @return State of the isothermal button
	 */
	public boolean isIsothermal() {
		return thermalButtons.isSelected("Isothermal");
	}

	/**
	 * Set the Adiabatic button to its selected state.
	 */
	public void setAdiabatic() {
        thermalButtons.setSelected("Adiabatic");
	}

	/**
	 * @return State of the adiabatic button
	 */
	public boolean isAdiabatic() {
        return thermalButtons.isSelected("Adiabatic");
	}

	/**
	 * Sets an action to be performed after a isothermal/adiabatic button is
	 * pressed and the integrator's isothermality has been set.
	 */
	public void setRadioGroupPostAction(final IAction action) {
	    thermalButtons.setPostAction(new IAction() {
	        public void actionPerformed() {
	            configureSliderAccessibility();
	            action.actionPerformed();
	        }
	    });
	}

	/**
	 * Set the current value for the temperature slider/text box.
	 */
    public void setTemperature(double value) {
        temperatureSlider.setValue(value);
    }

	/**
	 * @return  Current value of the temperature slider/text box.
	 */
	public double getTemperature() {
		return temperatureSlider.getValue();
	}

	/**
	 * Set whether the temperature text box should be displayed.
	 */
    public void setShowValues(boolean b) {
    	temperatureSlider.setShowValues(b);
    }

	/**
	 * Set whether the temperature text box should be editable.
	 */
    public void setEditValues(boolean b) {
    	temperatureSlider.setEditValues(b);
    }

	/**
	 * Set the minimum value for the temperature.
	 */
    public void setMinimum(double min) {
         temperatureSlider.setMinimum(min);
    }

	/**
	 * Set the maximum value for the temperature.
	 */
    public void setMaximum(double max) {
        temperatureSlider.setMaximum(max);
    }

	/**
	 * Set the number of "major" values that should be shown on the
	 * temperature slider.
	 */
    public void setSliderMajorValues(int major) {
    	temperatureSlider.setNMajor(major);
    }

    public void setSliderMinorValues(int minor) {
        temperatureSlider.setNMinor(minor);
    }

    /**
     * @return The panel that holds all graphical objects for the DeviceThermoSlider.
     */
    public Component graphic(Object obj) {
    	return temperaturePanel;
    }

	/**
	 * Set the unit of temperature.
	 */
    public void setUnit(Unit u) {
    	temperatureSlider.setUnit(u);
        String suffix = (u.symbol().length() > 0) ? " ("+u.symbol()+")" : "";
        temperaturePanel.setBorder(new TitledBorder(null, "Set Temperature" + suffix, TitledBorder.CENTER, TitledBorder.TOP));
    }

    /**
     * Set the precision of the scrollbar.
     * @param prec Number of significant digits after the "dot".
     */
    public void setPrecision(int prec) {
    	temperatureSlider.setPrecision(prec);
    }

	/**
	 * Set the temperature modifier object.
	 */
    public void setModifier(Modifier mod) {
        temperatureSlider.setModifier(mod);
    }

    /**
     * @return Temperature value modifier.
     */
    public Modifier getModifier() {
        return temperatureSlider.getModifier();
        
    }

	/**
	 * Set the temperature slider controller.
	 */
    public void setController(Controller cont) {
    	super.setController(cont);
        temperatureSlider.setController(cont);
        thermalButtons.setController(cont);
    }

	/**
	 * Set the post slider value changed action.
	 */
    public void setSliderPostAction(IAction action) {
        temperatureSlider.setPostAction(action);
    }

    //
    //main method to test device
    //
    public static void main(String[] args) {
        final String APP_NAME = "Device Thermo Slider";

        
        etomica.space.Space sp = etomica.space3d.Space3D.getInstance();
        etomica.simulation.Simulation sim = new etomica.simulation.Simulation(sp);
        final SimulationGraphic graphic = new SimulationGraphic(sim, APP_NAME, sp, sim.getController());

        DeviceThermoSlider device = new DeviceThermoSlider(new Controller(), new IntegratorVelocityVerlet(null, sim.getRandom(), 1, 1, sp));
        device.setMinimum(100.0);
        device.setMaximum(1000.0);
        device.setTemperature(250.0);

        graphic.getPanel().controlPanel.remove(graphic.getController().graphic());
        graphic.add(device);
        graphic.makeAndDisplayFrame(APP_NAME);

    }

    private void configureSliderAccessibility() {
        if (isAdiabatic()) {
        	temperatureSlider.getSlider().setEnabled(false);
        	temperatureSlider.getTextField().setEnabled(false);
        }
        else {
        	temperatureSlider.getSlider().setEnabled(true);
        	temperatureSlider.getTextField().setEnabled(true);
        }		
	}
}
