/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.vle;

import etomica.action.IAction;
import etomica.action.activity.Controller;
import etomica.graphics.Device;
import etomica.graphics.DeviceSlider;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorManagerMC;
import etomica.modifier.Modifier;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.units.Unit;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Temperature;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * Custom temperature slider for GEMC simulations.  The only real difference is
 * that this connects to the 2 sub-integrators (GEMC integrator is not an
 * IntegratorBox).
 * 
 * @author Andrew Schultz
 */
public class DeviceThermoSliderGEMC extends Device {

	private JPanel        temperaturePanel;  // main panel for thermo device PRIVATE
	private DeviceSlider  temperatureSlider; // Do not make make accessible
	private JRadioButton  buttonAdiabatic;   // Do not make make accessible
	private JRadioButton  buttonIsothermal;  // Do not make make accessible
	protected IntegratorBox    integrator1, integrator2;
    protected IntegratorManagerMC integrator3;

	private final int DEFAULT_MIN_TEMPERATURE = 0;
	private final int DEFAULT_MAX_TEMPERATURE = 300;

	public DeviceThermoSliderGEMC(Controller cont) {

        //adiabatic/isothermal radio button
        ButtonGroup thermalGroup = new ButtonGroup();
        buttonAdiabatic = new JRadioButton("Adiabatic");
        buttonIsothermal = new JRadioButton("Isothermal");
        thermalGroup.add(buttonAdiabatic);
        thermalGroup.add(buttonIsothermal);

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

        // Tie the isothermal/adiabatic setting to the selectable status of
        // temperature slider
        ToggleButtonListener myListener = new ToggleButtonListener();
        buttonAdiabatic.setSelected(true);
        buttonAdiabatic.addActionListener(myListener);
        buttonIsothermal.addActionListener(myListener);

        temperaturePanel = new JPanel(new GridBagLayout());
        temperaturePanel.setBorder(new TitledBorder(null, "Set Temperature", TitledBorder.CENTER, TitledBorder.TOP));
        GridBagConstraints gbc1 = new GridBagConstraints();
        gbc1.gridx = 0;  gbc1.gridy = 1;
        gbc1.gridwidth = 1;
        temperaturePanel.add(buttonAdiabatic, gbc1);
        gbc1.gridx = 1;  gbc1.gridy = 1;
        gbc1.gridwidth = 1;
        temperaturePanel.add(buttonIsothermal,gbc1);
        gbc1.gridx = 0;  gbc1.gridy = 2;
        gbc1.gridwidth = 2;
        temperaturePanel.add(temperatureSlider.graphic(),gbc1);
    }

	public void setIsothermalButtonsVisibility(boolean doShowIsothermalButtons) {
	    buttonIsothermal.setVisible(doShowIsothermalButtons);
        buttonAdiabatic.setVisible(doShowIsothermalButtons);
	}

	public boolean getIsothermalButtonsVisibility() {
	    return buttonIsothermal.isVisible();
	}
	
	/**
	 * Set the Isothermal button to its selected state.
	 */
	public void setIsothermal() {
		buttonIsothermal.setSelected(true);
		radioButtonChangeByClient();
		configureSliderAccessibility();
	}

	/**
	 * @return State of the isothermal button
	 */
	public boolean isIsothermal() {
		return buttonIsothermal.isSelected();
	}

	/**
	 * Set the Adiabatic button to its selected state.
	 */
	public void setAdiabatic() {
		buttonAdiabatic.setSelected(true);
		radioButtonChangeByClient();
		configureSliderAccessibility();
	}

	/**
	 * @return State of the adiabatic button
	 */
	public boolean isAdiabatic() {
		return buttonAdiabatic.isSelected();
	}

	private void radioButtonChangeByClient() {
		if(integrator1 != null && integrator2 != null) {
	        controller.doActionNow(integratorBoxIsoChangeSetIso);
	    }
	}

	/**
	 * Add the specified listener to the list of listeners that
	 * will get invoked when the temperature slider  value changes.
	 * @param listener
	 */
	public void addTemperatureSliderListener(ChangeListener listener) {
		temperatureSlider.getSlider().addChangeListener(listener);
	}

	/**
	 * Add the specified listener to the list of listeners that
	 * will get invoked when the isothermal or adiabatic radio button
	 * is pushed.
	 * @param listener
	 */
	public void addRadioGroupActionListener(ActionListener listener) {
		buttonAdiabatic.addActionListener(listener);
		buttonIsothermal.addActionListener(listener);
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

    public void doUpdate() {
        temperatureSlider.doUpdate();
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
        if (integrator1 != null && integrator2 != null) {
            // invoke setIntegartor again so that the isothermal/adiabatic
            // listener gets updated
            setIntegrators(integrator1, integrator2, integrator3);
        }
    }

	/**
	 * Set the post slider value changed action.
	 */
    public void setSliderPostAction(IAction action) {
        temperatureSlider.setPostAction(action);
    }

    /**
     * Sets the integrator for the device.  Adds actions to the device's
     * controller to inform the integrator when the temperature and
     * isothermal/adiabatic selection has changed based upon the type
     * of integrator passed in.
     * @param newIntegrator1
     * @param newIntegrator2
     * @param newIntegrator3
     */
    public void setIntegrators(IntegratorBox newIntegrator1, IntegratorBox newIntegrator2, IntegratorManagerMC newIntegrator3) {
        integrator1 = newIntegrator1;
    	integrator2 = newIntegrator2;
        integrator3 = newIntegrator3;
        temperatureSlider.setModifier(new Modifier() {

            public Dimension getDimension() {
                return Temperature.DIMENSION;
            }

            public String getLabel() {
                return "Temperature";
            }

            public double getValue() {
                return integrator1.getTemperature();
            }

            public void setValue(double newValue) {
                integrator1.setTemperature(newValue);
                integrator2.setTemperature(newValue);
                integrator3.setTemperature(newValue);
            }
    	});

    	ActionListener actionListen = new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
				controller.doActionNow(integratorBoxIsoChangeSetIso);
            }
        };

        addRadioGroupActionListener(actionListen);
        if (integrator1.isIsothermal()) {
            setIsothermal();
        }
        else {
            setAdiabatic();
        }
    }

    //
    //main method to test device
    //
    public static void main(String[] args) {
        final String APP_NAME = "Device Thermo Slider";

        Controller controller = new Controller();
        DeviceThermoSliderGEMC device = new DeviceThermoSliderGEMC(controller);
        device.setMinimum(100.0);
        device.setMaximum(1000.0);
        device.setTemperature(250.0);

        Space sp = Space2D.getInstance();
        etomica.simulation.Simulation sim = new etomica.simulation.Simulation(sp);
        final SimulationGraphic graphic = new SimulationGraphic(sim, APP_NAME);
        graphic.getPanel().controlPanel.remove(graphic.getController().graphic());
        graphic.add(device);
        graphic.makeAndDisplayFrame(APP_NAME);

    }

    private IAction integratorBoxIsoChangeSetIso = new IAction() {
        public void actionPerformed() {
            integrator1.setIsothermal(isIsothermal());
            integrator2.setIsothermal(isIsothermal());
        }
    };

    private void configureSliderAccessibility() {
        if(buttonAdiabatic.isSelected()) {
        	temperatureSlider.getSlider().setEnabled(false);
        	temperatureSlider.getTextField().setEnabled(false);
        }
        else {
        	temperatureSlider.getSlider().setEnabled(true);
        	temperatureSlider.getTextField().setEnabled(true);
        }		
	}

    /**
     * Private class that toggles the state of the temperature slider and
     * temperature text box based on the adiabatic/isothermal button currently
     * selected.  The slider/text box is selectable under isothermal conditions
     * and unselectable when adiabatic is selected.
     *
     */
    private class ToggleButtonListener implements ActionListener {
    	public void actionPerformed(ActionEvent e) {
    		configureSliderAccessibility();
        }
    }

}
