/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.action.IAction;
import etomica.action.controller.Controller;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.modifier.ModifyAction;
import etomica.units.Unit;
import etomica.units.systems.UnitSystem;
import etomica.util.Strings;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;

/**
 * Device the changes a property using a graphical spinner, via a Modifier.
 *
 * @author Andrew Schultz
 * @see ModifierGeneral
 */
public class DeviceSpinner extends Device {
    
    /**
     * Descriptive text label to be displayed with the value
     * No longer used -- instead apply label via a title border on the spinner itself
     */
    private String label;
    /**
     * Modifier connecting the spinner to the property
     */
    protected ModifyAction modifyAction;
    /**
     * Actual spinner model
     */
    protected SpinnerNumberModel spinnerModel;
    /**
     * Actual JSpinner
     */
    protected JSpinner spinner;
    /**
     * Object with property being modulated
     */
    protected Object component;
    /**
     * Property being modulated
     */
    protected String property;
    /**
     * Values for maximum and minimum values of spinner in double form
     */
    private double minimum, maximum;

    private JPanel panel; 
    
    /** 
     * horizontal position of TitledBorder.  Default is TitledBorder.LEFT
     */    
    private int borderAlignment = TitledBorder.LEFT;
    
    /**
     * Layout instance to show spinner 
     */    
    private GridBagLayout gbLayout;    
    private GridBagConstraints gbConst; 
    private boolean showBorder = false;
    protected IAction targetAction;

    public DeviceSpinner(Controller controller) {
        super(controller);
        init();
    }
    
    /**
     * Constructs a spinner connected to the given property of the given object
     */
    public DeviceSpinner(Controller controller, Object object, String property) {
        this(controller, new ModifierGeneral(object, property));
        component = object;
        this.property = property;
    }
    /**
     * Constructs a spinner connected to the get/set Value methods of the given Modifier
     */
    public DeviceSpinner(Controller controller, Modifier m) {
        this(controller);
        //set component and property in some way
        setModifier(m);
    }

    private void init() {
        spinnerModel = new SpinnerNumberModel(0, 0, 500, 1);
        spinner = new JSpinner(spinnerModel);
        spinner.setFont(new java.awt.Font("",0,15));
        gbLayout = new GridBagLayout();
        gbConst = new GridBagConstraints();    
        panel = new JPanel();  
        setLabel("");
        panel.setLayout(gbLayout);
        setMinimum(0);
        setMaximum(500);
        spinner.addChangeListener(new SpinnerListener());

        gbConst.gridx = 0; gbConst.gridy = 0;
        gbLayout.setConstraints(spinner, gbConst);
        panel.add(spinner);
    }
    
    
    /**
     * Override superclass setUnit method to update label when unit is changed
     */
    public void setUnit(Unit u) {
        super.setUnit(u);
        setLabelDefault();
    }

    public void setBorderAlignment(int align) {
    	borderAlignment = align;
    }

    public void setModifier(Modifier m) {
        if(m == null) throw new NullPointerException();
        modifyAction = null;
        if(unit == null) {
            setUnit(m.getDimension().getUnit(UnitSystem.SIM));
        }
        spinnerModel.setValue((int)unit.fromSim(m.getValue()));
        modifyAction = new ModifyAction(m);
        //targetActions needs to be distinct from modifyAction, in case subclasses want to do more than modifyAction when slider is moved
        targetAction = modifyAction;
        setLabelDefault();
        setMinimum(getMinimum());
        setMaximum(getMaximum());
    }
    public final Modifier getModifier() {return modifyAction.getWrappedModifier();}
    
    public void doUpdate() {
        double value = unit.fromSim(modifyAction.getValue());
        spinnerModel.setValue((int)value);
    }
    
    public String getProperty() {return property;}
    public void setProperty(String s) {
        property = s;
        if(component != null) setModifier(new ModifierGeneral(component, property));
    }
    
    public Object getComponent() {return component;}
    public void setComponent(Object obj) {
        component = obj;
        if(property != null) setModifier(new ModifierGeneral(component,property));
    }
    
    public double getMinimum() {return minimum;}
    /**
     * Sets minimum value of spinner
     */
    public void setMinimum(double min) {
        minimum = min;
        ModifyAction tmpModifier = modifyAction;
        modifyAction = null;
        spinnerModel.setMinimum((int)min);
        modifyAction = tmpModifier;
    }
        
    public double getMaximum() {return maximum;}
    /**
     * Sets maximum value of slider
     */
    public void setMaximum(double max) {
        maximum = max;
        ModifyAction tmpModifier = modifyAction;
        modifyAction = null;
        spinnerModel.setMaximum((int)max);
        modifyAction = tmpModifier;
    }
    
    public double getValue(){return spinnerModel.getNumber().doubleValue();}    
    public void setValue(double d){
        //updating the spinner will update the text field
        spinnerModel.setValue((int)d);
    }
    
    /**
     * Returns the GUI element for display in the simulation.
     */
    public java.awt.Component graphic() {
        return panel;
    }

    /**
     * Sets the value of a descriptive label using the meter's label and the unit's symbol (abbreviation).
     */
    private void setLabelDefault() {
        String suffix = (unit.symbol().length() > 0) ? " ("+unit.symbol()+")" : "";
        if(modifyAction != null) 
            setLabel(Strings.capitalize(modifyAction.getWrappedModifier().getLabel())+suffix);
    }

    /**
     * Sets the value of a descriptive label using the given string.
     */
    public void setLabel(String s){
        label = s;
        if(s == null || s.equals("") || !showBorder) panel.setBorder(new javax.swing.border.EmptyBorder(0,0,0,0));
        else {
        	TitledBorder border = new TitledBorder(s);
        	border.setTitleJustification(borderAlignment);
        	panel.setBorder(border);
        }
    }
    
    /**
     * @return the current instance of the descriptive label.
     */
    public String getLabel() {return label;}
    
    public void setShowBorder(boolean b) {
        showBorder = b;
        setLabel(label);
    }
    public boolean isShowBorder() {return showBorder;}
    
    /**
     * @return a handle to the SpinnerNumberModel instance used by this slider device
     */
    public SpinnerNumberModel getSpinner() {return spinnerModel;}
    
    /**
     * @return a handle to the JPanel instance used by this slider device
     */    
    public JPanel getPanel() { return panel; }    
     
     /**
      * Spinner listener, which relays the spinner change events to the modifier
      */
     protected class SpinnerListener implements ChangeListener, java.io.Serializable {
         public void stateChanged(ChangeEvent evt) {
             double newValue = spinnerModel.getNumber().intValue();
             if(modifyAction!=null) {
                 modifyAction.setValueForAction(unit.toSim(newValue));
                 try {
                     doAction(targetAction);
                 }
                 catch (RuntimeException e) {
                     //XXX We want to put the slider back where it was (which is hopefully
                     //what the backend element has now), but changing the 
                     //slider position fires an event that brings us back here.  
                     //null-out modifyAction to prevent recursion!
                     ModifyAction actualModifyAction = modifyAction;
                     modifyAction = null;
                     // new value is the old value!
                     newValue = unit.fromSim(actualModifyAction.getValue());
                     spinnerModel.setValue((int)newValue);
                     modifyAction = actualModifyAction;
                 }
             }
         }
     }
}
