package simulate;
import simulate.gui.SimEditorTabMenu;
import java.beans.PropertyEditorSupport;
import java.util.LinkedList;

/**
 * Editor to create a list of instantiated objects of type MeterAbstract.
 */
public class MeterEditor extends PropertyEditorSupport implements java.io.Serializable {
    
    Object[] meterObjects;
    String[] meterNames;
    
    public MeterEditor() {
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
        int nMeters = SimEditorTabMenu.getMeterEditor().getComponentList().size()+1;
        meterNames = new String[nMeters];
        meterObjects = new Object[nMeters];
        meterNames[0] = "";
        meterObjects[0] = null;
        for(int i=1; i<nMeters; i++) {
            meterObjects[i] = ((MeterAbstract)SimEditorTabMenu.getMeterEditor().getComponentList().elementAt(i-1));
            meterNames[i] = meterObjects[i].toString();
        }
    }
}