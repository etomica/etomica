package etomica;
import etomica.units.*;
import java.beans.PropertyEditor;
import java.beans.PropertyEditorSupport;
import javax.swing.JComboBox;
import javax.swing.JTextField;

/**
 * Editor for a property of type Unit.
 */
public class UnitEditor extends PropertyEditorSupport 
            implements java.awt.event.ItemListener, 
                java.io.Serializable,
                javax.swing.JComboBox.KeySelectionManager {

    
    private Unit unit;
    private Prefix prefix;
    private BaseUnit baseUnit;
    private Class[] availableUnits;
    private StringBuffer[] availableUnitNames;
    private JComboBox unitList;
    
    public UnitEditor(Unit currentUnit) {
        super();
        setupUnits(currentUnit);
    }
    
    private void setupUnits(Unit currentUnit) {
        Dimension dimension;
        if(currentUnit == null) dimension = Dimension.NULL;
        if(currentUnit instanceof UnitRatio) dimension = Dimension.NULL; //chokes on UnitRatio because cannot instantiate in setupNames method
        else dimension = currentUnit.dimension();
        unit = currentUnit;
        prefix = unit.prefix();
        baseUnit = unit.baseUnit();
        availableUnits = BaseUnit.all(dimension); //all BaseUnit classes of the given dimensions
        availableUnitNames = new StringBuffer[availableUnits.length];
        for(int i=0; i<availableUnitNames.length; i++) {
            availableUnitNames[i] = new StringBuffer();
        }
        setupNames();
        if(unitList != null) {
            unitList.removeItemListener(this);
            unitList.removeAllItems();
            for(int i=0; i<availableUnitNames.length; i++) {
                unitList.addItem(availableUnitNames[i]);
            }
        }
        else {
            unitList = new JComboBox(availableUnitNames);
        }
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
    
    public Unit getUnit() {return unit;}
    
    public void setValue(Object obj) {
        super.setValue(obj);
        setupUnits((Unit)obj);
    }
    
    /**
     * ItemListener interface implementation.
     */
     public void itemStateChanged(java.awt.event.ItemEvent evt) {
        Class baseClass = availableUnits[unitList.getSelectedIndex()];
	    try {
	        baseUnit = (BaseUnit)baseClass.newInstance();
	    }
	    catch(InstantiationException e) {System.out.println(e.toString()); System.exit(1);}
	    catch(IllegalAccessException e) {System.out.println(e.toString()); System.exit(1);}
	    unit = new Unit(prefix, baseUnit);
	    setValue(unit);
	 }
	 
	 /**
	  * Implementation of KeySelectionManager interface that uses key press to
	  * select a unit prefix (instead of selecting a combo box item).
	  */
	 public int selectionForKey(char aKey, javax.swing.ComboBoxModel aModel) {
	    Prefix newPrefix = Prefix.keySelect(aKey);
	    if(newPrefix != null) { //update unit if an appropriate key was pressed
	        prefix = newPrefix;
	        unit = new Unit(prefix, baseUnit);
    	    setupNames();
	        setValue(unit);
	        unitList.getTopLevelAncestor().repaint();
	    }
	    return -1;  //indicate no selection of combobox items
	 }
	    
	 /**
	  * Returns a combo box that can be used to select the unit.
	  */
	 public JComboBox unitSelector() {return unitList;}
	 
    
    /**
     * Returns value in simulation units.
     */
    public String getJavaInitializationString() {
        return getValue().toString();
    }
}