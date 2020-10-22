/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

      package etomica.normalmode;

import etomica.action.IAction;
import etomica.action.controller.Controller;
import etomica.graphics.Device;
import etomica.graphics.DeviceSlider;
import etomica.graphics.SimulationGraphic;
import etomica.modifier.Modifier;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Null;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;


/**
 * @author taitan
 */
public class DeviceEigenvaluesSlider extends Device {

	private JPanel eValNumPanel;  // main panel for eigenvalues # device PRIVATE
	private DeviceSlider eValNumSlider; // Do not make accessible
	private JRadioButton  buttonAllEVal;   // Do not make accessible
	private JRadioButton  buttonOneEVal;  // Do not make accessible
	private final IntegratorHarmonic integrator;

	private final int DEFAULT_MIN_nEIGENVECTORS = 0;
	private final int DEFAULT_MAX_nEIGENVECTORS = 24;

	public DeviceEigenvaluesSlider(Controller cont, IntegratorHarmonic integrator) {
		super(cont);

		this.integrator = integrator;
		eValNumSlider = new DeviceSlider(controller);
		eValNumSlider.setModifier(new Modifier(){

			public Dimension getDimension() {
				return Null.DIMENSION;
			}

			public String getLabel() {
				return "Eigenvalues #";
			}

			public double getValue() {
				return integrator.getEValNum();
			}

			public void setValue(double eValValue) {
				integrator.setEValNum((int)eValValue);

			}

		});

		ActionListener actionListen = evt -> controller.submitActionInterrupt(integratorBoxChangeSetOneEval);

		addRadioGroupActionListener(actionListen);

		if (this.integrator.isOneEVal()) {
			setOneEVal();

		}
		else {
			setAllEVal();
		}


		//using all eigenvalues or individual radio button
		ButtonGroup eValGroup = new ButtonGroup();
		buttonAllEVal = new JRadioButton("All Normal Modes");
        buttonOneEVal = new JRadioButton("One Normal Mode");
        eValGroup.add(buttonAllEVal);
        eValGroup.add(buttonOneEVal);

        //eigenvalues selector
        eValNumSlider.setShowValues(true);
        eValNumSlider.setEditValues(true);
        eValNumSlider.setMinimum(DEFAULT_MIN_nEIGENVECTORS);
        eValNumSlider.setMaximum(DEFAULT_MAX_nEIGENVECTORS);
        eValNumSlider.setNMajor(6);
        eValNumSlider.setValue(0);
        eValNumSlider.getSlider().setEnabled(false);
        eValNumSlider.getTextField().setEnabled(false);

        // Tie the "all eigenvalues"/"one eigenvalue" setting to the selectable status of
        // eigenvalues slider
        ToggleButtonListener myListener = new ToggleButtonListener();
        buttonAllEVal.setSelected(true);
        buttonAllEVal.addActionListener(myListener);
        buttonOneEVal.addActionListener(myListener);

        eValNumPanel = new JPanel(new GridBagLayout());
        eValNumPanel.setBorder(new TitledBorder(null, "Set Normal Modes", TitledBorder.CENTER, TitledBorder.TOP));
        GridBagConstraints gbc1 = new GridBagConstraints();
        gbc1.gridx = 0;  gbc1.gridy = 1;
        gbc1.gridwidth = 1;
        eValNumPanel.add(buttonAllEVal, gbc1);
        gbc1.gridx = 0;  gbc1.gridy = 2;
        gbc1.gridwidth = 1;
        eValNumPanel.add(buttonOneEVal,gbc1);
        gbc1.gridx = 0;  gbc1.gridy = 3;
        gbc1.gridwidth = 2;
        eValNumPanel.add(eValNumSlider.graphic(),gbc1);
        
    }

	public void setOneEValButtonsVisibility(boolean doShowOneEValButtons) {
	    buttonOneEVal.setVisible(doShowOneEValButtons);
        buttonAllEVal.setVisible(doShowOneEValButtons);
	}

	public boolean getOneEValButtonsVisibility() {
	    return buttonOneEVal.isVisible();
	}
	
	/**
	 * Set the "eigenvalues" button to its selected state.
	 */
	public void setOneEVal() {
		buttonOneEVal.setSelected(true);
		radioButtonChangeByClient();
		configureSliderAccessibility();
	}

	/**
	 * @return State of the 'eigenvalues' button
	 */
	public boolean isOneEVal() {
		return buttonOneEVal.isSelected();
	}

	/**
	 * Set the "All eigenvalues" button to its selected state.
	 */
	public void setAllEVal() {
		buttonAllEVal.setSelected(true);
		radioButtonChangeByClient();
		configureSliderAccessibility();
	}

	/**
	 * @return State of the "All Eigenvalues" button
	 */
	public boolean isAllEVal() {
		return buttonAllEVal.isSelected();
	}

	private void radioButtonChangeByClient() {
		if(integrator != null) {
	        controller.submitActionInterrupt(integratorBoxChangeSetOneEval);
	    }
	}

	/**
	 * Add the specified listener to the list of listeners that
	 * will get invoked when the 'eigenvalues #' slider value changes.
	 * @param listener
	 */
	public void addEValNumSliderListener(ChangeListener listener) {
		eValNumSlider.getSlider().addChangeListener(listener);
	}

	/**
	 * Add the specified listener to the list of listeners that
	 * will get invoked when the "Eigenvalue" or "All Eigenvalues" radio button
	 * is pushed.
	 * @param listener
	 */
	public void addRadioGroupActionListener(ActionListener listener) {
		buttonAllEVal.addActionListener(listener);
		buttonOneEVal.addActionListener(listener);
	}

	/**
	 * Set the current value for the eigenvalue # slider/text box.
	 */
	
    public void setEValNum(int value) {
        eValNumSlider.setValue(value);
    }

	/**
	 * @return  Current value of the eigenvalue # slider/text box.
	 */
	public double getEValNum() {
		return eValNumSlider.getValue();
	}

	/**
	 * Set whether the eigenvalue # text box should be displayed.
	 */
    public void setShowEValvalues(boolean b) {
    	eValNumSlider.setShowValues(b);
    }

	/**
	 * Set whether the eigenvalue # text box should be editable.
	 */
    public void setEditValues(boolean b) {
    	eValNumSlider.setEditValues(b);
    }

	/**
	 * Set the minimum value for the eigenvalue #.
	 */
    public void setMinimum(int min) {
         eValNumSlider.setMinimum(min);
    }

	/**
	 * Set the maximum value for the eigenvalue #.
	 * 
	 *  max-1 (because the value in the array start from 0)
	 */
    public void setMaximum(int max) {
    	eValNumSlider.setMaximum(max-1);
    }

	/**
	 * Set the number of "major" values that should be shown on the
	 * eigenvalue slider.
	 */
    public void setSliderMajorValues(int major) {
    	eValNumSlider.setNMajor(major);
    }

    /**
	 * @return The panel that holds all graphical objects for the DeviceEigenvaluesSlider.
	 */
	public Component graphic() {
		return eValNumPanel;
	}


	/**
	 * Set the eigenvalue # modifier object.
	 */
    public void setModifier(Modifier mod) {
        eValNumSlider.setModifier(mod);
    }

    /**
     * @return eigenvalue # value modifier.
     */
    public Modifier getModifier() {
        return eValNumSlider.getModifier();
        
    }

	/**
	 * Set the post slider value changed action.
	 */
    public void setSliderPostAction(IAction action) {
    	eValNumSlider.setPostAction(action);
    }


    private final IAction integratorBoxChangeSetOneEval = new IAction() {
        public void actionPerformed() {
        	integrator.setOneEVal(isOneEVal());
        	
        }
    };

    private void configureSliderAccessibility() {
        if(buttonAllEVal.isSelected()) {
        	eValNumSlider.getSlider().setEnabled(false);
        	eValNumSlider.getTextField().setEnabled(false);
        }
        else {
        	eValNumSlider.getSlider().setEnabled(true);
        	eValNumSlider.getTextField().setEnabled(true);
        }		
	}

    /**
     * Private class that toggles the state of the eigenvalue # slider and
     * temperature text box based on the "All eigenvalues"/"One eigenvalue" button currently
     * selected.  The slider/text box is selectable under "One eigenvalue" conditions
     * and unselectable when "All eigenvalue" is selected.
     *
     */
    private class ToggleButtonListener implements ActionListener {
    	public void actionPerformed(ActionEvent e) {
    		configureSliderAccessibility();
        }
    }
    
}
