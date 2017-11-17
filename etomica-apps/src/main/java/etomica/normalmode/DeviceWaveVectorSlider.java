/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

      package etomica.normalmode;

import java.awt.Component;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.ButtonGroup;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeListener;

import etomica.action.IAction;
import etomica.action.activity.Controller;
import etomica.graphics.Device;
import etomica.graphics.DeviceSlider;
import etomica.graphics.SimulationGraphic;
import etomica.modifier.Modifier;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Null;


/**
 * @author taitan
 *
 */
public class DeviceWaveVectorSlider extends Device {

	private JPanel        waveVectorNumPanel;  // main panel for wave vectors device PRIVATE
	private DeviceSlider  waveVectorNumSlider; // Do not make make accessible
	private JRadioButton  buttonAllWV;   // Do not make make accessible
	private JRadioButton  buttonOneWV;  // Do not make make accessible
	protected IntegratorHarmonic integrator;

	private final int DEFAULT_MIN_nWAVEVECTORS = 0;
	private final int DEFAULT_MAX_nWAVEVECTORS = 100;

	public DeviceWaveVectorSlider(Controller cont) {
		
        //using all wave vectors or individual radio button
        ButtonGroup waveVectorGroup = new ButtonGroup();
        buttonAllWV = new JRadioButton("All Wave Vectors");
        buttonOneWV = new JRadioButton("Wave Vector");
        waveVectorGroup.add(buttonAllWV);
        waveVectorGroup.add(buttonOneWV);

        //waveVector selector
        waveVectorNumSlider = new DeviceSlider(controller);
        waveVectorNumSlider.setShowValues(true);
        waveVectorNumSlider.setEditValues(true);
        waveVectorNumSlider.setMinimum(DEFAULT_MIN_nWAVEVECTORS);
        waveVectorNumSlider.setMaximum(DEFAULT_MAX_nWAVEVECTORS);
        waveVectorNumSlider.setNMajor(4);
        waveVectorNumSlider.setValue(0);
        waveVectorNumSlider.getSlider().setEnabled(false);
        waveVectorNumSlider.getTextField().setEnabled(false);

        setController(cont);

        // Tie the "all wave vectors"/"one wave vector" setting to the selectable status of
        // waveVector slider
        ToggleButtonListener myListener = new ToggleButtonListener();
        buttonAllWV.setSelected(true);
        buttonAllWV.addActionListener(myListener);
        buttonOneWV.addActionListener(myListener);

        waveVectorNumPanel = new JPanel(new GridBagLayout());
        waveVectorNumPanel.setBorder(new TitledBorder(null, "Set Wave Vector", TitledBorder.CENTER, TitledBorder.TOP));
        GridBagConstraints gbc1 = new GridBagConstraints();
        gbc1.gridx = 0;  gbc1.gridy = 1;
        gbc1.gridwidth = 1;
        waveVectorNumPanel.add(buttonAllWV, gbc1);
        gbc1.gridx = 1;  gbc1.gridy = 1;
        gbc1.gridwidth = 1;
        waveVectorNumPanel.add(buttonOneWV,gbc1);
        gbc1.gridx = 0;  gbc1.gridy = 2;
        gbc1.gridwidth = 2;
        waveVectorNumPanel.add(waveVectorNumSlider.graphic(),gbc1);
    }

	public void setOneWVButtonsVisibility(boolean doShowOneWVButtons) {
	    buttonOneWV.setVisible(doShowOneWVButtons);
        buttonAllWV.setVisible(doShowOneWVButtons);
	}

	public boolean getOneWVButtonsVisibility() {
	    return buttonOneWV.isVisible();
	}
	
	/**
	 * Set the "wave vector" button to its selected state.
	 */
	public void setOneWV() {
		buttonOneWV.setSelected(true);
		radioButtonChangeByClient();
		configureSliderAccessibility();
	}

	/**
	 * @return State of the 'wave vector' button
	 */
	public boolean isOneWV() {
		return buttonOneWV.isSelected();
	}

	/**
	 * Set the "All Wave Vectors" button to its selected state.
	 */
	public void setAllWV() {
		buttonAllWV.setSelected(true);
		radioButtonChangeByClient();
		configureSliderAccessibility();
	}

	/**
	 * @return State of the "All Wave Vectors" button
	 */
	public boolean isAllWV() {
		return buttonAllWV.isSelected();
	}

	private void radioButtonChangeByClient() {
		if(integrator != null) {
	        controller.doActionNow(integratorBoxChangeSetOneWV);
	    }
	}

	/**
	 * Add the specified listener to the list of listeners that
	 * will get invoked when the 'wave vectors #' slider value changes.
	 * @param listener
	 */
	public void addWaveVectorNumSliderListener(ChangeListener listener) {
		waveVectorNumSlider.getSlider().addChangeListener(listener);
	}

	/**
	 * Add the specified listener to the list of listeners that
	 * will get invoked when the "Wave Vector" or "All Wave Vectors" radio button
	 * is pushed.
	 * @param listener
	 */
	public void addRadioGroupActionListener(ActionListener listener) {
		buttonAllWV.addActionListener(listener);
		buttonOneWV.addActionListener(listener);
	}

	/**
	 * Set the current value for the wavevector # slider/text box.
	 */
	
    public void setWaveVectorNum(int value) {
        waveVectorNumSlider.setValue(value);
    }

	/**
	 * @return  Current value of the wavevector # slider/text box.
	 */
	public double getWaveVectorNum() {
		return waveVectorNumSlider.getValue();
	}

	/**
	 * Set whether the wavevector # text box should be displayed.
	 */
    public void setShowValues(boolean b) {
    	waveVectorNumSlider.setShowValues(b);
    }

	/**
	 * Set whether the wavevector # text box should be editable.
	 */
    public void setEditValues(boolean b) {
    	waveVectorNumSlider.setEditValues(b);
    }

	/**
	 * Set the minimum value for the wavevector #.
	 */
    public void setMinimum(int min) {
         waveVectorNumSlider.setMinimum(min);
    }

	/**
	 * Set the maximum value for the wavevector #.
	 * 
	 *  max-1 (because the value in the array start from 0)
	 */
    public void setMaximum(int max) {
    	waveVectorNumSlider.setMaximum(max-1);
    }

	/**
	 * Set the number of "major" values that should be shown on the
	 * temperature slider.
	 */
    public void setSliderMajorValues(int major) {
    	waveVectorNumSlider.setNMajor(major);
    }

    /**
     * @return The panel that holds all graphical objects for the DeviceThermoSlider.
     */
    public Component graphic(Object obj) {
    	return waveVectorNumPanel;
    }


	/**
	 * Set the wavevector # modifier object.
	 */
    public void setModifier(Modifier mod) {
        waveVectorNumSlider.setModifier(mod);
    }

    /**
     * @return wavevector # value modifier.
     */
    public Modifier getModifier() {
        return waveVectorNumSlider.getModifier();
        
    }

	/**
	 * Set the wavevector # slider controller.
	 */
   
    public void setController(Controller cont) {
    	super.setController(cont);
    	waveVectorNumSlider.setController(cont);
        
    	if (integrator != null) {
            setIntegrator(integrator);
        }
    }

	/**
	 * Set the post slider value changed action.
	 */
    public void setSliderPostAction(IAction action) {
    	waveVectorNumSlider.setPostAction(action);
    }

    /**
     * Sets the integrator for the device.  Adds actions to the device's
     * controller to inform the integrator when the wavevector # and
     * "All Wave Vectors"/"One Wave Vector" selection has changed based upon the type
     * of integrator passed in.
     * @param i Integrator
     */
    
    public void setIntegrator(IntegratorHarmonic i) {
    	integrator = i;
    	waveVectorNumSlider.setModifier(new Modifier(){

			public Dimension getDimension() {
				return Null.DIMENSION;
			}

			public String getLabel() {
				return "Wave Vector #";
			}

			public double getValue() {
				return integrator.getWaveVectorNum();
			}

			public void setValue(double newValue) {
				integrator.setWaveVectorNum((int)newValue);
				
			}
    		
    	});

    	ActionListener actionListen = new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
				controller.doActionNow(integratorBoxChangeSetOneWV);
            }
        };

        addRadioGroupActionListener(actionListen);
        
        if (i.isOneWV()) {
        	setOneWV();
        	
        }
        else {
            setAllWV();
        }
        
    }


    
    //
    //main method to test device
    //
    public static void main(String[] args) {
        final String APP_NAME = "Device Wave Vectors Number Slider";

       
        etomica.space.Space sp = etomica.space1d.Space1D.getInstance();
        etomica.simulation.Simulation sim = new etomica.simulation.Simulation(sp);
        
        DeviceWaveVectorSlider device = new DeviceWaveVectorSlider(new Controller());
        device.setMinimum(0);
        device.setMaximum(100);
        device.setWaveVectorNum(0);
        
        
        final SimulationGraphic graphic = new SimulationGraphic(sim, APP_NAME, sp, sim.getController());
        graphic.getPanel().controlPanel.remove(graphic.getController().graphic());
        graphic.add(device);
        graphic.makeAndDisplayFrame(APP_NAME);

    }
    
    
    
    
    
    
    
    

    
    
    private IAction integratorBoxChangeSetOneWV = new IAction() {
        public void actionPerformed() {
        	integrator.setOneWV(isOneWV());
        	
        }
    };

    private void configureSliderAccessibility() {
        if(buttonAllWV.isSelected()) {
        	waveVectorNumSlider.getSlider().setEnabled(false);
        	waveVectorNumSlider.getTextField().setEnabled(false);
        }
        else {
        	waveVectorNumSlider.getSlider().setEnabled(true);
        	waveVectorNumSlider.getTextField().setEnabled(true);
        }		
	}

    /**
     * Private class that toggles the state of the wave vectors # slider and
     * temperature text box based on the "All Wave Vectors"/"One Wave Vector" button currently
     * selected.  The slider/text box is selectable under "One Wave Vector" conditions
     * and unselectable when "All Wave Vectors" is selected.
     *
     */
    private class ToggleButtonListener implements ActionListener {
    	public void actionPerformed(ActionEvent e) {
    		configureSliderAccessibility();
        }
    }



}
