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
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;


/**
 * @author taitan
 *
 */
public class DeviceCellNumXYSlider extends Device {

	private JPanel        numCellPanel;  // main panel for cell number device PRIVATE
	private DeviceSlider  xCellNumSlider; 
	private DeviceSlider yCellNumSlider; // Do not make make accessible
	private JRadioButton  buttonXComp;   // Do not make make accessible
	private JRadioButton  buttonYComp;  // Do not make make accessible
	
	private final int DEFAULT_MIN_nCells = 1;
	private final int DEFAULT_MAX_nCells = 50;

    protected ISpecies species;
    protected Box box;
    
	
	public DeviceCellNumXYSlider(Controller cont) {
		
        //using x-axis or y-axis radio button
        ButtonGroup numCellGroup = new ButtonGroup();
        buttonXComp = new JRadioButton("x-Cell");
        buttonYComp = new JRadioButton("y-Cell");
        numCellGroup.add(buttonXComp);
        numCellGroup.add(buttonYComp);

        //x-CellNum selector
        xCellNumSlider = new DeviceSlider(controller);
        xCellNumSlider.setShowValues(true);
        xCellNumSlider.setEditValues(true);
        xCellNumSlider.setMinimum(DEFAULT_MIN_nCells);
        xCellNumSlider.setMaximum(DEFAULT_MAX_nCells);
        xCellNumSlider.setNMajor(4);
        xCellNumSlider.setValue(0);
        xCellNumSlider.getSlider().setEnabled(true);
        xCellNumSlider.getTextField().setEnabled(true);

        //y-CellNum selector
        yCellNumSlider = new DeviceSlider(controller);
        yCellNumSlider.setShowValues(true);
        yCellNumSlider.setEditValues(true);
        yCellNumSlider.setMinimum(DEFAULT_MIN_nCells);
        yCellNumSlider.setMaximum(DEFAULT_MAX_nCells);
        yCellNumSlider.setNMajor(4);
        yCellNumSlider.setValue(0);
        yCellNumSlider.getSlider().setEnabled(false);
        yCellNumSlider.getTextField().setEnabled(false);
        
        setController(cont);

        // Tie the "x-axis"/"y-axis" setting to the selectable status of
        // numCells slider
        ToggleButtonListener myListener = new ToggleButtonListener();
        buttonXComp.setSelected(true);
        buttonXComp.addActionListener(myListener);
        buttonYComp.addActionListener(myListener);

        numCellPanel = new JPanel(new GridBagLayout());
        numCellPanel.setBorder(new TitledBorder(null, "Set x- and y- Cell Numbers", TitledBorder.CENTER, TitledBorder.TOP));
        GridBagConstraints gbc1 = new GridBagConstraints();
        
        gbc1.gridx = 0;  gbc1.gridy = 1;
        gbc1.gridwidth = 1;
        numCellPanel.add(buttonXComp, gbc1);
        gbc1.gridx = 1;  gbc1.gridy = 1;
        gbc1.gridwidth = 1;
        numCellPanel.add(buttonYComp,gbc1);
        
        gbc1.gridx = 0;  gbc1.gridy = 2;
        gbc1.gridwidth = 2;
        numCellPanel.add(xCellNumSlider.graphic(),gbc1);
        gbc1.gridx = 0;  gbc1.gridy = 3;
        gbc1.gridwidth = 3;
        numCellPanel.add(yCellNumSlider.graphic(),gbc1);
    }
	
	/**
	 * Set the "y-axis" button to its selected state.
	 */
	public void setYComp() {
		buttonYComp.setSelected(true);
		configureSliderAccessibility();
	}

	/**
	 * @return State of the 'y-axis' button
	 */
	public boolean isYComp() {
		return buttonYComp.isSelected();
	}

	/**
	 * Set the "x-axis" button to its selected state.
	 */
	public void setXComp() {
		buttonXComp.setSelected(true);
		configureSliderAccessibility();
	}

	/**
	 * @return State of the "x-axis" button
	 */
	public boolean isXComp() {
		return buttonXComp.isSelected();
	}

	
	
	/**
	 * Add the specified listener to the list of listeners that
	 * will get invoked when the 'xCell #' and 'yCell #' slider value changes.
	 * @param listener
	 */
	public void addXYCellNumSliderListener(ChangeListener listener) {
		xCellNumSlider.getSlider().addChangeListener(listener);
		yCellNumSlider.getSlider().addChangeListener(listener);
	}

	/**
	 * Add the specified listener to the list of listeners that
	 * will get invoked when the "x-axis" or "y-axis" radio button
	 * is pushed.
	 * @param listener
	 */
	public void addRadioGroupActionListener(ActionListener listener) {
		buttonXComp.addActionListener(listener);
		buttonYComp.addActionListener(listener);
	}

	
	/////////////////////////////////////////////////////////////////////////////
	/**
	 * Set the current value for the xCell # slider/text box.
	 */
	
    public void setXCellNum(int value) {
        xCellNumSlider.setValue(value);
    }

	/**
	 * @return  Current value of the xCell # slider/text box.
	 */
	public double getXCellNum() {
		return xCellNumSlider.getValue();
	}
	
	/////////////////////////////////////////////////////////////////////////////
	/**
	 * Set the current value for the yCell # slider/text box.
	 */
	
    public void setYCellNum(int value) {
        yCellNumSlider.setValue(value);
    }

	/**
	 * @return  Current value of the yCell # slider/text box.
	 */
	public double getYCellNum() {
		return yCellNumSlider.getValue();
	}
	/////////////////////////////////////////////
	

	/**
	 * Set whether the x- and y-Cell # text box should be displayed.
	 */
    public void setShowValues(boolean b) {
    	xCellNumSlider.setShowValues(b);
    	yCellNumSlider.setShowValues(b);
    }

	/**
	 * Set whether the x- and y-Cell # text box should be editable.
	 */
    public void setEditValues(boolean b) {
    	xCellNumSlider.setEditValues(b);
    	yCellNumSlider.setEditValues(b);
    }

	/**
	 * Set the minimum value for the x- and y-Cell #.
	 */
    public void setMinimum(int min) {
         xCellNumSlider.setMinimum(min);
         yCellNumSlider.setMinimum(min);
    }

	/**
	 * Set the maximum value for the x- and y-Cell #.
	 * 
	 */
    public void setMaximum(int max) {
    	xCellNumSlider.setMaximum(max);
    	yCellNumSlider.setMaximum(max/2);
    }

	/**
	 * Set the number of "major" values that should be shown on the
	 * x- and y- CellNum slider.
	 */
    public void setSliderMajorValues(int major) {
    	xCellNumSlider.setNMajor(major);
    	yCellNumSlider.setNMajor(major/2);
    }

    /**
     * @return The panel that holds all graphical objects for the DeviceCellNumXYSlider.
     */
    public Component graphic() {
        return numCellPanel;
    }


	/**
	 * Set the x-Cell # modifier object.
	 */
    public void setXCellModifier(Modifier mod) {
        xCellNumSlider.setModifier(mod);
    }

    /**
     * @return x-Cell # value modifier.
     */
    public Modifier getXCellModifier() {
        return xCellNumSlider.getModifier();
        
    }

    
	/**
	 * Set the y-Cell # modifier object.
	 */
    public void setYCellModifier(Modifier mod) {
        yCellNumSlider.setModifier(mod);
    }

    /**
     * @return y-Cell # value modifier.
     */
    public Modifier getYCellModifier() {
        return yCellNumSlider.getModifier();
        
    }
    
	/**
	 * Set the x- and y- Cell # slider controller.
	 */
    public void setController(Controller cont) {
    	super.setController(cont);
    	xCellNumSlider.setController(cont);
        yCellNumSlider.setController(cont);
       
    }

	/**
	 * Set the post slider value changed action.
	 */
    public void setXSliderPostAction(IAction action) {
       	setYCellModifier(new ModifierYCells2D(box, species, (int)getXCellNum()));
    	xCellNumSlider.setPostAction(action);
    	setXCellModifier(new ModifierXCells2D(box, species, (int)getYCellNum()));
    }

    public void setYSliderPostAction(IAction action) {
    	setXCellModifier(new ModifierXCells2D(box, species, (int)getYCellNum()));
    	yCellNumSlider.setPostAction(action);
    	setYCellModifier(new ModifierYCells2D(box, species, (int)getXCellNum()));
    }
    
    public void setBox(Box newBox) {
        box = newBox;
        if (species != null) {
            init();
        }
    }
    
    public void setSpecies(ISpecies newSpecies) {
        species = newSpecies;
        if (box != null) {
            init();
        }
    }
    
    public Box getBox() {
        return box;
    }
    
    public ISpecies getSpecies() {
        return species;
    }

    
    protected void init(){
    	
    	setXCellModifier(new ModifierXCells2D(box, species, (int)getYCellNum()));
       	setYCellModifier(new ModifierYCells2D(box, species, (int)getXCellNum()));
    	
    }

    private void configureSliderAccessibility() {
        if(buttonXComp.isSelected()) {
        	xCellNumSlider.getSlider().setEnabled(true);
        	xCellNumSlider.getTextField().setEnabled(true);
        	
        	yCellNumSlider.getSlider().setEnabled(false);
        	yCellNumSlider.getTextField().setEnabled(false);
        	init();
        }
        else {
        	yCellNumSlider.getSlider().setEnabled(true);
        	yCellNumSlider.getTextField().setEnabled(true);
        	
        	xCellNumSlider.getSlider().setEnabled(false);
        	xCellNumSlider.getTextField().setEnabled(false);
        	init();
        }		
	}

    /**
     * Private class that toggles the state of the cell # slider and
     * cell # text box based on the "x-axis"/"y-axis" button currently
     * selected.  The y- slider/text box is selectable under "y-axis" conditions
     * and x- slider/text box is selectable when "x-axis" is selected.
     *
     */
    private class ToggleButtonListener implements ActionListener {
    	public void actionPerformed(ActionEvent e) {
    		configureSliderAccessibility();
        }
    }




    

}
