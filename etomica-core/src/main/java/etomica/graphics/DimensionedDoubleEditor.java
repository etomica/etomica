/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;
import java.beans.PropertyEditor;
import java.beans.PropertyEditorSupport;

import javax.swing.JComboBox;

import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Null;
import etomica.units.Prefix;
import etomica.units.PrefixedUnit;
import etomica.units.SimpleUnit;
import etomica.units.Unit;
import etomica.units.systems.UnitSystem;

/**
 * Editor to set the value of a double-type property having units associated with it.
 * Holds an instance of a Unit class that is applied to convert the value
 * to and from simulation units when setting and getting the value.  Different
 * forms of the set and get methods interpret input/output value in either 
 * simulation units, or the units associated with this editor, as indicated
 * by the comments for each method.
 */
 
 //might want to revise to incorporate UnitEditor, which was built
 //using the units-handling features of this editor
public class DimensionedDoubleEditor extends PropertyEditorSupport 
            implements java.awt.event.ItemListener, 
                java.io.Serializable,
                javax.swing.JComboBox.KeySelectionManager {

    // holds value in editor units, and permits display/editing as a string
    private PropertyEditor valueEditor;
    
//    private JTextField editor;
    private PrefixedUnit unit;
    private Prefix prefix;
    private Unit baseUnit;
    private Class[] availableUnits;
    private StringBuffer[] availableUnitNames;
    private JComboBox unitList;
    
    public DimensionedDoubleEditor(Dimension dimension) {
        super();
        valueEditor = java.beans.PropertyEditorManager.findEditor(Double.TYPE);
        if(dimension == null) dimension = Null.DIMENSION;
        Unit defaultUnit = dimension.getUnit(UnitSystem.SIM);
        unit = (defaultUnit instanceof PrefixedUnit) ? (PrefixedUnit)defaultUnit : new PrefixedUnit(defaultUnit);
        prefix = unit.prefix();
        baseUnit = unit.unit();
        //this is broken  availableUnits = SimpleUnit.all(dimension); //all BaseUnit classes of the given dimensions
        availableUnitNames = new StringBuffer[availableUnits.length];
        for(int i=0; i<availableUnitNames.length; i++) {
            availableUnitNames[i] = new StringBuffer();
        }
        setupNames();
        unitList = new JComboBox(availableUnitNames);
        //set unitlist combo box so it shows current value of unit
        //do this before adding listeners
        for(int i=0; i<availableUnits.length; i++) {
            if(availableUnitNames[i].toString().equals(unit.toString())) {
                unitList.setSelectedIndex(i);
            }
        }
        //add unitList listeners
        unitList.addItemListener(this);
        unitList.setKeySelectionManager(this);
    }
    
    private void setupNames() {
        SimpleUnit base = null;
        for(int i=0; i<availableUnits.length; i++) {
	        try {
	            base = (SimpleUnit)availableUnits[i].newInstance();
	        }
	        catch(InstantiationException e) {System.out.println(e.toString()); System.exit(1);}
	        catch(IllegalAccessException e) {System.out.println(e.toString()); System.exit(1);}
            if(base != null) {
                availableUnitNames[i].replace(0,availableUnitNames[i].length(),
                                            new PrefixedUnit(prefix, base).toString());
            }
            else {
                availableUnitNames[i] = null;
            }
        }
    }
    
    public Unit getUnit() {return unit;}
    
    /**
     * ItemListener interface implementation.
     */
     public void itemStateChanged(java.awt.event.ItemEvent evt) {
        Object value = getValue();
        Class baseClass = availableUnits[unitList.getSelectedIndex()];
	    try {
	        baseUnit = (SimpleUnit)baseClass.newInstance();
	    }
	    catch(InstantiationException e) {System.out.println(e.toString()); System.exit(1);}
	    catch(IllegalAccessException e) {System.out.println(e.toString()); System.exit(1);}
	    unit = new PrefixedUnit(prefix, baseUnit);
	    setValue(value);
	 }
	 
	 /**
	  * Implementation of KeySelectionManager interface that uses key press to
	  * select a unit prefix (instead of selecting a combo box item).
	  */
	 public int selectionForKey(char aKey, javax.swing.ComboBoxModel aModel) {
        Object value = getValue();
        String sKey;
        if(aKey == ' ') sKey = "";
        else if(aKey == 'D') sKey = "da";
        else sKey = Character.toString(aKey);
	    Prefix newPrefix = Prefix.keySelect(sKey);
	    if(newPrefix != null) { //update unit if an appropriate key was pressed
	        prefix = newPrefix;
	        unit = new PrefixedUnit(prefix, baseUnit);
    	    setupNames();
	        setValue(value);
	        unitList.getTopLevelAncestor().repaint();
	    }
	    return -1;  //indicate no selection of combobox items
	 }
	    
	 
	 /**
	  * Returns a combo box that can be used to select the unit.
	  */
	 public JComboBox unitSelector() {return unitList;}
	 
     public PropertyEditor valueEditor() {
//        if(editor==null) editor = new etomica.gui.PropertyText(this);
        return valueEditor;
     }   

    /**
     * Returns the value in the editor's units, formatted as a text string.
     */
    public String getAsText() {
        return valueEditor.getAsText();
    }
    
    /**
     * Sets the value, interpreted as being in the editor's units, from a text string.
     */
    public void setAsText(String s) {
        valueEditor.setAsText(s);
        firePropertyChange();
    }
    
    /**
     * Returns the value in simulation units.
     */
    public Object getValue() {
        double value = ((Double)valueEditor.getValue()).doubleValue();  //value in editor units
        value = unit.toSim(value);  //convert to simulation units
        return new Double(value);   //return value in simulation units
    }
    
    /**
     * Sets the value to the given value, which is taken to be in simulation units.
     */
    public void setValue(Object obj) {
          //value in simulation units
        double value = (obj!=null) ? ((Double)obj).doubleValue() : Double.NaN;
        value = unit.fromSim(value);                 //convert to editor units
        valueEditor.setValue(new Double(value));     //store in editor units
        firePropertyChange();
    }
    
    /**
     * Sets the value to the given value, which is taken to be in simulation units.
     * Does not fire a property-change event.  This is useful for displays that want
     * to constantly update the current value without causing all the follow-up usually
     * performed in response to the change.  Particularly useful with updating meter
     * value displays, where the update comes from within the simulation thread, not 
     * by user input.
     */
    public void setValue(double value) {
        valueEditor.setValue(new Double(unit.fromSim(value)));
    }
    
    /**
     * Returns value in simulation units.
     */
    public String getJavaInitializationString() {
    //    String value = getAsText();
    //   String unitName = unit.getClass().getName();
    //    return unitName+".toSim("+value+")";
        return getValue().toString();
    }
}
