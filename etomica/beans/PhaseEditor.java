package simulate;
import simulate.*;
import simulate.gui.SimEditorTabMenu;
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
        int nPhases = SimEditorTabMenu.getPhaseEditor().getComponentList().size();
        phaseNames = new String[nPhases];
        phaseObjects = new Object[nPhases];
        for(int i=0; i<SimEditorTabMenu.getPhaseEditor().getComponentList().size(); i++) {
            phaseObjects[i] = ((Phase)SimEditorTabMenu.getPhaseEditor().getComponentList().elementAt(i));
            phaseNames[i] = phaseObjects[i].toString();
        }
    }
}
