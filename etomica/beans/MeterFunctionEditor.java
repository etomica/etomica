package simulate;
import simulate.gui.SimEditorTabMenu;
import java.beans.PropertyEditorSupport;

/**
 * Editor to create a list of instantiated objects of type MeterFunction.
 */
public class MeterFunctionEditor extends PropertyEditorSupport implements java.io.Serializable {
    
    Object[] meterObjects;
    String[] meterNames;
    
    public MeterFunctionEditor() {
        super();
    }

    public String getAsText() {
        makeMeterLists();
        if(meterNames.length == 0 || getValue() == null) return "null";
        return getValue().toString();
    }
    
    public void setAsText(String s) {
        makeMeterLists();
        for(int i=0; i<meterNames.length; i++) {
            if(meterNames[i].equals(s)) {
                setValue(meterObjects[i]);
                firePropertyChange();
                return;
            }
        }
    }
    
    public String[] getTags() {
        makeMeterLists();
        return meterNames;
    }

    protected void makeMeterLists() {
        javax.swing.DefaultListModel list = SimEditorTabMenu.getMeterEditor().getComponentList();
        //count how many of all MeterAbstract objects are MeterFunction objects
        int nMeters = 1;
        for(int i=0; i<list.size(); i++) {
            if(list.elementAt(i) instanceof MeterFunction) nMeters++;
        }
        //loop again, save all the MeterFunction objects in the arrays
        meterNames = new String[nMeters];
        meterObjects = new Object[nMeters];
        meterNames[0] = "";
        meterObjects[0] = null;
        nMeters = 1;
        for(int i=0; i<list.size(); i++) {
            if(list.elementAt(i) instanceof MeterFunction) {
                meterObjects[nMeters] = list.elementAt(i);
                meterNames[nMeters] = meterObjects[nMeters].toString();
                nMeters++;
            }
        }
    }
}