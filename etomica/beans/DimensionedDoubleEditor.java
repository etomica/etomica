package simulate;
import simulate.units.*;
import simulate.gui.SimEditorTabMenu;
import java.beans.PropertyEditor;
import java.beans.PropertyEditorSupport;
import javax.swing.JComboBox;

/**
 * Editor to set the value of a double-type property having units associated with it.
 * Holds an instance of a Unit class that is applied to convert the value
 * to and from simulation units when setting and getting the value.  Different
 * forms of the set and get methods interpret input/output value in either 
 * simulation units, or the units associated with this editor, as indicated
 * by the comments for each method.
 */
public class DimensionedDoubleEditor extends PropertyEditorSupport 
            implements java.awt.event.ItemListener, 
                java.io.Serializable,
                javax.swing.JComboBox.KeySelectionManager {
    
    private PropertyEditor valueEditor;
    private Unit unit;
    private Prefix prefix;
    private BaseUnit baseUnit;
    private Class[] availableUnits;
    private StringBuffer[] availableUnitNames;
    private JComboBox unitList;
    
    public DimensionedDoubleEditor(Dimension dimension) {
        super();
        valueEditor = java.beans.PropertyEditorManager.findEditor(Double.TYPE);
        unit = dimension.defaultIOUnit();
        prefix = unit.prefix();
        baseUnit = unit.baseUnit();
        availableUnits = BaseUnit.all(dimension); //all BaseUnit classes of the given dimensions
        availableUnitNames = new StringBuffer[availableUnits.length];
        for(int i=0; i<availableUnitNames.length; i++) {
            availableUnitNames[i] = new StringBuffer();
        }
        setUpNames();
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
    
    private void setUpNames() {
        BaseUnit base = null;
        for(int i=0; i<availableUnits.length; i++) {
	        try {
	            base = (BaseUnit)availableUnits[i].newInstance();
	        }
	        catch(InstantiationException e) {System.out.println(e.toString()); System.exit(1);}
	        catch(IllegalAccessException e) {System.out.println(e.toString()); System.exit(1);}
            if(base != null) {
                availableUnitNames[i].replace(0,availableUnitNames[i].length(),
                                            new Unit(prefix, base).toString());
            }
            else {
                availableUnitNames[i] = null;
            }
        }
    }
    
    /**
     * ItemListener interface implementation.
     */
     public void itemStateChanged(java.awt.event.ItemEvent evt) {
        Object value = getValue();
        Class baseClass = availableUnits[unitList.getSelectedIndex()];
	    try {
	        baseUnit = (BaseUnit)baseClass.newInstance();
	    }
	    catch(InstantiationException e) {System.out.println(e.toString()); System.exit(1);}
	    catch(IllegalAccessException e) {System.out.println(e.toString()); System.exit(1);}
	    unit = new Unit(prefix, baseUnit);
	    setValue(value);
	 }
	 
	 /**
	  * Implementation of KeySelectionManager interface that uses key press to
	  * select a unit prefix (instead of selecting a combo box item).
	  */
	 public int selectionForKey(char aKey, javax.swing.ComboBoxModel aModel) {
        Object value = getValue();
	    Prefix newPrefix = Prefix.keySelect(aKey);
	    if(newPrefix != null) { //update unit if an appropriate key was pressed
	        prefix = newPrefix;
	        unit = new Unit(prefix, baseUnit);
    	    setUpNames();
	        setValue(value);
	    }
	    return -1;  //indicate no selection of combobox items
	 }
	    
	 
	 /**
	  * Returns a combo box that can be used to select the unit.
	  */
	 public JComboBox unitSelector() {return unitList;}
        

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
}