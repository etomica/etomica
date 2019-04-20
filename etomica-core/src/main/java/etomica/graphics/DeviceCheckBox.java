/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

//This class includes a main method to demonstrate its use
package etomica.graphics;

import etomica.modifier.ModifierBoolean;
import etomica.modifier.ModifyBooleanAction;

import javax.swing.*;

/**
 * Button that toggles a boolean value via a check box.  This device can connect to any object
 * capable of switching between two states.  The device operates through a 
 * ModifierBoolean instance that must be connected to the state of the
 * controlled object.
 * 
 * @author David Kofke
 */
public class DeviceCheckBox extends Device {

    protected ModifyBooleanAction modifyAction;
    private ModifierBoolean modifier;
    private JCheckBox box;
    private boolean currentValue = false;
    
    public DeviceCheckBox(ModifierBoolean modifier) {
        this("Select", modifier);
    }
    public DeviceCheckBox(String label, ModifierBoolean modifier) {
        super();
        init(label, modifier);
    }
    
    private void init(String label, ModifierBoolean modifier) {
        this.modifier = modifier;
        if (modifier != null) {
            modifyAction = new ModifyBooleanAction(modifier);
            currentValue = modifier.getBoolean();
        }
        box = new JCheckBox(label,currentValue);
        box.addActionListener( new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                toggle();
            }
        });
    }

    /**
     * Sets the button to the given state.
     */
    public void setState(boolean b) {
        if(b != currentValue) toggle();
    }
    /**
     * Returns the currents true/false state of the button.
     */
    public boolean getState() {return currentValue;}
    
    /**
     * Toggles the state of the button (true to false, false to true).
     */
    public void toggle() {
        currentValue = !currentValue;
        modifyAction.setValueForAction(currentValue);
        doAction(modifyAction);
        box.setSelected(currentValue);
    }
    
    /**
     * Specifies the boolean modifier that is set between true and false by the button.
     */
    public ModifierBoolean getModifier() {return modifier;}
    /**
     * Returns the boolean modifier set by this button.
     */
    public void setModifier(ModifierBoolean newModifier) {
        modifier = newModifier;
        modifier.setBoolean(currentValue);
        modifyAction = new ModifyBooleanAction(modifier);
        currentValue = modifier.getBoolean();
    }
    
    /**
     * Returns the GUI box element for display in the simulation.
     */
    public java.awt.Component graphic(Object obj) {
        return box;
    }
    
    
    /**
     * @return a handle to the JCheckBox instance used by this device
     */
    public JCheckBox getCheckBox() {return box;}
}
