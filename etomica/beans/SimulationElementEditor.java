package etomica;
import java.beans.*;
//Java2 imports
//import java.util.LinkedList;
//import java.util.Iterator;

import etomica.utility.LinkedList;
import etomica.utility.Iterator;

/**
 * A property editor for setting the a Phase property in an object.
 */

public class SimulationElementEditor extends PropertyEditorSupport {
    
    Object[] elementObjects;
    String[] elementNames;
    Class elementClass;
    
    public SimulationElementEditor(Class elementClass) {
        super();
        this.elementClass = elementClass;
        makeElementLists();
    }

    public String getAsText() {
        makeElementLists();
        if(elementNames.length == 0 || getValue() == null) return "null";
        return getValue().toString();
    }
    
    public void setAsText(String s) {
        makeElementLists();
        for(int i=0; i<elementNames.length; i++) {
            if(elementNames[i].equals(s)) {
                setValue(elementObjects[i]);
                firePropertyChange();
                return;
            }
        }
    }

    public String[] getTags() {return elementNames;}

    protected void makeElementLists() {
        LinkedList elementList = etomica.Simulation.instance.elementList(elementClass);
        elementNames = new String[elementList.size()];
        elementObjects = new Object[elementList.size()];
        int i=0;
        for(Iterator iter=elementList.iterator(); iter.hasNext(); i++) {
            elementObjects[i] = iter.next();
            elementNames[i] = elementObjects[i].toString();
        }
    }
}
