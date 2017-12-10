/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;
import java.awt.Color;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.KeyEvent;
import java.util.Enumeration;
import java.util.Vector;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import etomica.modifier.Modifier;
import etomica.modifier.ModifyAction;
import etomica.units.dimensions.Null;
import etomica.units.systems.UnitSystem;
import etomica.util.Constants;
import etomica.util.EnumeratedType;

/**
 * A simple device the permits editing of a single value via a textbox 
 * with an associated label.
 *
 * @author David Kofke
 */
 
public class DeviceBox extends Device implements javax.swing.event.ChangeListener, DeviceBoxValueChangedListener {
    
    /**
     * Descriptive text label to be displayed with the value
     */
    protected JLabel label;
    protected String labelString;
    private Constants.CompassDirection labelPosition = Constants.CompassDirection.NORTH;
    /**
     * Object for displaying the value as a text field
     */
    protected JTextField textField;
    /**
     * Displayed panel that holds the label and value
     * (not yet used; meant to implement to make lightweight display)
     */
    protected JPanel panel;
    /**
     * Modifier connecting the slider to the property
     */
    protected ModifyAction modifyAction;

    /**
     * Integer specifying the number of significant figures to be displayed.
     * Default is 4.
     */
    protected int precision;

	private final java.awt.Insets borderInset = new java.awt.Insets(4, 4, 4, 4);

    private LabelType labelType;
    private boolean integer = false;
    private BoxListener listener = null;

    private Vector valueChangedListeners = new Vector();

    public DeviceBox() {
        super();
        panel = new JPanel(new java.awt.BorderLayout());
        label = null;
        labelString = null;
        textField = new JTextField("");
        textField.setEditable(true);
        panel.add(textField, java.awt.BorderLayout.CENTER);
        setLabelType(LabelType.STRING);
        unit = Null.UNIT;
        setPrecision(4);
    }

    /** 
     * calls doUpdate method.  Implementation of ChangeListener interface.
     */
    public void stateChanged(javax.swing.event.ChangeEvent evt) {
        doUpdate();
    }
 
    /**
     * Updates the display of the box with the current value given by the modifier.
     */
    public void doUpdate() {
        if(modifyAction == null) {
        	textField.setText("");
        	return;
        }
        if(integer) {
            textField.setText(Integer.toString((int)unit.fromSim(modifyAction.getWrappedModifier().getValue())));
        }
        else {
            textField.setText(DisplayTextBox.format(unit.fromSim(modifyAction.getWrappedModifier().getValue()),precision));
        }
    }
    
    public void setEditable(boolean b) {textField.setEditable(b);}
    public boolean isEditable() {return textField.isEditable();}
    
    /**
     * Accessor method to set the physical units of the displayed value.
     * Text describing unit is used in label.
     */
    public void setUnit(etomica.units.Unit u) {
        unit = u;
        doUpdate();
        setLabel(labelString);
    }
    
    public java.awt.Component graphic(Object obj) {return panel;}
    
    /**
     * Accessor method of the precision, which specifies the number of significant figures to be displayed.
     */
    public int getPrecision() {return precision;}

    /**
     * Accessor method of the precision, which specifies the number of significant figures to be displayed.
     */
    public void setPrecision(int n) {
        textField.setColumns(n+2);
        precision = n;
    }
    
    /**
     * Sets a flag indicating if the value should be displayed as an integer.
     */
    public void setInteger(boolean b) {integer = b;}
    public boolean isInteger() {return integer;}
    
    /**
     * Specifies the modifier that receives the edit.
     */
    public void setModifier(Modifier m) {
        if(m == null) {
        	modifyAction = null;

        	// Don't want to listen for events in textField if
        	// there is not modifier.
        	if(listener != null) {
                textField.removeActionListener(listener);
                textField.removeFocusListener(listener);
                textField.removeKeyListener(listener);
        	}

        	// If there is no modifier, make the textField uneditable
        	textField.setEditable(false);

        	listener = null;
        	throw new NullPointerException();
        }
        if(unit == null) {
            setUnit(m.getDimension().getUnit(UnitSystem.SIM));
        }
        modifyAction = new ModifyAction(m);
        setLabel(labelString);
        doUpdate();

        listener = new BoxListener();

        textField.addActionListener(listener);
        textField.addKeyListener(listener);
        textField.addFocusListener(listener);
    }
    
    /**
     * Accessor method for the modifier that receives the edit.
     */
    public Modifier getModifier() {
        return modifyAction == null ? null : modifyAction.getWrappedModifier();
    }
    
    /**
     * Sets the value of a descriptive label using the given string.
     */
    public void setLabel(String s) {
        labelString = s;
        if(labelType == LabelType.BORDER) {
            panel.setBorder(new javax.swing.border.TitledBorder(getLabel()));
        }
        else if(labelType == LabelType.STRING) {
            label.setText(getLabel());
            setLabelPosition(labelPosition);
        }
    }
    
    /**
     * @return the current value of the descriptive label.
     */
    public String getLabel() {
        if (labelString != null) {
            return labelString;
        }
        if (modifyAction == null) {
            return "";
        }
        String symbolText = unit.symbol();
        return modifyAction.getWrappedModifier().getLabel()+
                ((symbolText.length() > 0) ? " ("+symbolText+")" : "");
    }
    
    /**
     * Sets the label type to "border" or "string".  With border, the label
     * text is part of the border.  With "string", the label is inside the 
     * border.
     */
    public void setLabelType(LabelType newLabelType) {
        if (newLabelType == labelType) {
            return;
        }
        labelType = newLabelType;
        if(labelType == LabelType.BORDER) {
            panel.remove(label);
            label = null;
            panel.setBorder(new javax.swing.border.TitledBorder(getLabel()));
        }
        else if(labelType == LabelType.STRING) {
            label = new JLabel(getLabel());
            panel.setBorder(new javax.swing.border.MatteBorder(borderInset, panel.getBackground()));
            setLabelPosition(labelPosition);
        }
    }

    public LabelType getLabelType() {
        return labelType;
    }

    public void setLabelPosition(Constants.CompassDirection position) {
        labelPosition = position;
        if(labelType != LabelType.STRING) return;
        panel.remove(label);
        panel.add(label,position.toString());//toString() returns the corresponding BorderLayout constant
        panel.revalidate();
        panel.repaint();
    }

    public void setTextBackground(Color color) {
    	textField.setBackground(color);
    }

    public Color getTextBackground() {
    	return textField.getBackground();
    }

    public void setBorderBackground(Color color) {
    	javax.swing.border.Border b = panel.getBorder();
    	javax.swing.border.MatteBorder border =
    		new javax.swing.border.MatteBorder(borderInset, color);

    	// Only set the border color if the border is a MatteBorder
    	if(b.getClass().isInstance(border)) {
    		panel.setBorder(border);
    	}
    	else {
            border = null;
        }

    }

    public Constants.CompassDirection getLabelPosition() {
        return labelPosition;
    }

    private class BoxListener extends java.awt.event.KeyAdapter implements java.awt.event.ActionListener, FocusListener {
        public void actionPerformed(java.awt.event.ActionEvent evt) {
            updateValue();
        }
        public void keyReleased(java.awt.event.KeyEvent evt) {
            if(!isInteger()) return;
            int key = evt.getKeyCode();
            int step = 0;
            if(key == KeyEvent.VK_UP || key == KeyEvent.VK_RIGHT) step = 1;
            else if(key == KeyEvent.VK_DOWN || key == KeyEvent.VK_LEFT) step = -1;
            if(step == 0) return;
            int i = Integer.parseInt(textField.getText()) + step;
            textField.setText(Integer.toString(i));
            updateValue();
        }
       
        public void focusGained(FocusEvent evt) {
        }
        public void focusLost(FocusEvent evt) {
            // user tabbed or clicked out of the text field
            updateValue();
       }
       
       /**
        * Take the value form the textfield and apply it to the simulation via 
        * the modifyAction.  Update the textbox accordingly. 
        */
       protected void updateValue() {
           double oldX = modifyAction.getWrappedModifier().getValue();
           double newX;
           try { 
               newX = unit.toSim(Double.parseDouble(textField.getText()));
           } catch(NumberFormatException ex) {//if bad format, just restore original value
               doUpdate();
               return;
           }
           if (newX != oldX) {
               modifyAction.setValueForAction(newX);
               try {
                   doAction(modifyAction);
               }
               catch (RuntimeException e) {
                   // if the modifier throws, we'll revert to the old value in doUpdate
               }
           }
           doUpdate();
           notifyValueChangedListeners();
       }
    }

    /**
     * Add an object to be notified when the DeviceBox value to changes
     */
    public void addValueChangedListener(DeviceBoxValueChangedListener dbListener) {
    	if(!valueChangedListeners.contains(dbListener)) {
    		valueChangedListeners.add(dbListener);
    	}
    } // end addValueChangedListener

    /**
     * Remove an object that is listening for the DeviceBox value to change
     */
    public void removeValueChangedListener(DeviceBoxValueChangedListener dbListener) {
    	if(valueChangedListeners.contains(dbListener)) {
    		valueChangedListeners.remove(dbListener);
    	}
    } // end removeValueChangedListener

    /**
     * Removes all objects that are listening for a value changed event.
     */
    public void removeAllValueChangedListeners() {
    	valueChangedListeners.clear();
    }

    public void deviceBoxValueChanged(DeviceBoxValueChangedEvent ev) {

        this.doUpdate();
    } // end deviceBoxValueChanged

    /**
     *  Will notify all registered objects that the DeviceBox value has changed.
     *  Notification takes place at the tail end of the updateValue() method.
     */
    private void notifyValueChangedListeners() {
    	Vector copyOfListeners = (Vector)(valueChangedListeners.clone());

    	int value = 0; //<NEW VALUE HERE>
    	DeviceBoxValueChangedEvent ev = new DeviceBoxValueChangedEvent(this, value);
    	Enumeration vcElements = copyOfListeners.elements();
    	while(vcElements.hasMoreElements()) {
    		DeviceBoxValueChangedListener listener =
    			(DeviceBoxValueChangedListener)vcElements.nextElement();
    		listener.deviceBoxValueChanged(ev);
    	}

    } // end private void notifyValueChangedListeners

    /**
     * Typed constant used to indicate the type of label to be used with the display.
     */
	public enum LabelType {
	    BORDER,
        STRING
    }


}
