package etomica;
import java.beans.PropertyEditorSupport;

/**
 * Editor for typed constants, such as VERTICAL/HORIZONTAL.
 * Each typed constant instance is one of a defined set of values that can be adopted
 * by a field of that type (e.g. a field of type Alignment can be HORIZONTAL or VERTICAL).
 * Each instance inherits the choices() method, which returns the complete set of values
 * that are available.  This editor returns this set in its getTags method, and it examines
 * string-formatted changes against the set when changing the value of the field it is editing.
 * @see Constants
 */
public class TypedConstantEditor extends PropertyEditorSupport implements java.io.Serializable {
    
    public static String version() {return "01.03.05";}
    
    Constants.TypedConstant[] choices;
    String[] labels;
    
    public TypedConstantEditor() {
        super();
    }
    public TypedConstantEditor(Constants.TypedConstant[] c) {
        super();
        setChoices(c);
    }
    
    public void setValue(Object obj) {
        super.setValue(obj);
        if(choices == null && obj != null) setChoices(((Constants.TypedConstant)obj).choices());
    }
    
    private void setChoices(Constants.TypedConstant[] c) {
        choices = c;
        if(choices == null) {
            labels = null; 
            return;
        }
        labels = new String[choices.length];
        for(int i=0; i<choices.length; i++) {
            labels[i] = choices[i].toString();
        }
    }        

    public String getAsText() {
        if(choices == null || choices.length == 0 || getValue() == null) return "null";
        return getValue().toString();
    }
    
    public void setAsText(String s) {
        for(int i=0; i<choices.length; i++) {
            if(labels[i].equals(s)) {
                setValue(choices[i]);
                firePropertyChange();
                return;
            }
        }
    }
    
    /**
     * Method examined by property sheet to get selections for combo box.
     */
    public String[] getTags() {
        return labels;
    }

}//end of class