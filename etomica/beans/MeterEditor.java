package simulate;
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
        java.util.LinkedList meterList = Simulation.instance.meterList();
        int nMeters = meterList.size()+1;
        meterNames = new String[nMeters];
        meterObjects = new Object[nMeters];
        meterNames[0] = "";
        meterObjects[0] = null;
        int i=1;
        for(java.util.Iterator iter=meterList.iterator(); iter.hasNext(); i++) {
            meterObjects[i] = (MeterAbstract)iter.next();
            meterNames[i] = meterObjects[i].toString();
        }
    }
}