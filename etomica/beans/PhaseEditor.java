package etomica;
import java.beans.*;
import java.util.*;

/**
 * A property editor for setting the a Phase property in an object.
 */

public class PhaseEditor extends PropertyEditorSupport {
    
    Object[] phaseObjects;
    String[] phaseNames;
    
    public PhaseEditor() {
        super();
        makePhaseLists();
    }

    public String getAsText() {
        makePhaseLists();
        if(phaseNames.length == 0 || getValue() == null) return "null";
        return getValue().toString();
    }
    
    public void setAsText(String s) {
        makePhaseLists();
        for(int i=0; i<phaseNames.length; i++) {
            if(phaseNames[i].equals(s)) {
                setValue(phaseObjects[i]);
                firePropertyChange();
                return;
            }
        }
    }

    public String[] getTags() {return phaseNames;}

    protected void makePhaseLists() {
        java.util.LinkedList phaseList = etomica.Simulation.instance.phaseList();
        phaseNames = new String[phaseList.size()];
        phaseObjects = new Object[phaseList.size()];
        int i=0;
        for(java.util.Iterator iter=phaseList.iterator(); iter.hasNext(); i++) {
            phaseObjects[i] = (Phase)iter.next();
            phaseNames[i] = phaseObjects[i].toString();
        }
    }
}
