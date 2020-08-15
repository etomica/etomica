/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

      package etomica.normalmode;

import etomica.action.IAction;
import etomica.action.activity.Controller;
import etomica.box.Box;
import etomica.graphics.Device;
import etomica.graphics.DeviceSlider;
import etomica.graphics.SimulationGraphic;
import etomica.modifier.Modifier;
import etomica.species.ISpecies;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeListener;
import java.awt.*;


/**
 * @author taitan
 *
 */
public class DeviceCellNum3DSlider extends Device {

	private JPanel        numCellPanel;  // main panel for cell number device PRIVATE
	private DeviceSlider  nCellNumSlider; 
	
	private final int DEFAULT_MIN_nCells = 1;
	private final int DEFAULT_MAX_nCells = 100;

    protected ISpecies species;
    protected Box box;
    
	
	public DeviceCellNum3DSlider(Controller cont) {
		
        //n-CellNum selector
        nCellNumSlider = new DeviceSlider(controller);
        nCellNumSlider.setShowValues(true);
        nCellNumSlider.setEditValues(true);
        nCellNumSlider.setMinimum(DEFAULT_MIN_nCells);
        nCellNumSlider.setMaximum(DEFAULT_MAX_nCells);
        nCellNumSlider.setNMajor(5);
        nCellNumSlider.getSlider().setEnabled(true);
        nCellNumSlider.getTextField().setEnabled(true);

        setController(cont);

        numCellPanel = new JPanel(new GridBagLayout());
        numCellPanel.setBorder(new TitledBorder(null, "Set 3D n-Cell Numbers", TitledBorder.CENTER, TitledBorder.TOP));
        GridBagConstraints gbc1 = new GridBagConstraints();

        gbc1.gridx = 0;  gbc1.gridy = 1;
        gbc1.gridwidth = 1;
        numCellPanel.add(nCellNumSlider.graphic(),gbc1);
    }
	
	/**
	 * Add the specified listener to the list of listeners that
	 * will get invoked when the 'nCell #' slider value changes.
	 * @param listener
	 */
	public void addCellNumSliderListener(ChangeListener listener) {
		nCellNumSlider.getSlider().addChangeListener(listener);
	}
	
	/**
	 * Set the current value for the N-Cell # slider/text box.
	 */
	
    public void setNCellNum(int value) {
        nCellNumSlider.setValue(value);
    }

	/**
	 * @return  Current value of the N-Cell # slider/text box.
	 */
	public double getNCellNum() {
		return nCellNumSlider.getValue();
	}
	

	/**
	 * Set the minimum value for the N-Cell #.
	 */
    public void setMinimum(int min) {
         nCellNumSlider.setMinimum(min);
    }

	/**
	 * Set the maximum value for the N-Cell #.
	 * 
	 */
    public void setMaximum(int max) {
    	nCellNumSlider.setMaximum(max);
    }
   
	/**
	 * Set the number of "major" values that should be shown on the
	 *	n- CellNum slider.
	 */
    public void setSliderMajorValues(int major) {
    	nCellNumSlider.setNMajor(major);
    }

    /**
     * @return The panel that holds all graphical objects for the DeviceCellNumXYSlider.
     */
    public Component graphic() {
        return numCellPanel;
    }


	/**
	 * Set the n-Cell # modifier object.
	 */
    public void setNCellModifier(Modifier mod) {
        nCellNumSlider.setModifier(mod);
    }

    /**
     * @return -Cell # value modifier.
     */
    public Modifier getNCellModifier() {
        return nCellNumSlider.getModifier();
        
    }

    
	/**
	 * Set the n- Cell # slider controller.
	 */
    public void setController(Controller cont) {
    	super.setController(cont);
    	nCellNumSlider.setController(cont);
    }

	/**
	 * Set the post slider value changed action.
	 */
    public void setNSliderPostAction(IAction action) {
    	setNCellModifier(new ModifierCells3D(box, species));
    	nCellNumSlider.setPostAction(action);
    }

    
    public void setBox(Box newBox) {
        box = newBox;
        if (species != null) {
        	setNCellModifier(new ModifierCells3D(box, species));
        }
    }
    
    public void setSpecies(ISpecies newSpecies) {
        species = newSpecies;
        if (box != null) {
        	setNCellModifier(new ModifierCells3D(box, species));
        }
    }
    
    public Box getBox() {
        return box;
    }
    
    public ISpecies getSpecies() {
        return species;
    }



}
