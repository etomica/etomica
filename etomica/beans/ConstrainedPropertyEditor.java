package etomica;
import etomica.*;
import java.beans.*;
import java.awt.*;

/**
 * No-op editor for a constrained or otherwise uneditable property.
 * Holds the place of an editor to avoid have a null for the editor of such properties.
 *
 * @author David Kofke
 */
public class ConstrainedPropertyEditor extends PropertyEditorSupport implements java.io.Serializable {
    
    //Value of property for reference by editor, but without attempting to 
    //fire property change if altered
    private Object localValue;
    
    public String getAsText() {return null;}
    
    public boolean isPaintable() {return true;}
    
    public void paintValue(Graphics g, Rectangle box) { }
    
    public boolean supportsCustomEditor() {return false;}
    
    //Override so to not fire property change
    public Object getValue() {return localValue;}
    public void setValue(Object obj) {localValue = obj;}
    
}