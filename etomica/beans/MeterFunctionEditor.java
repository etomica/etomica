package etomica;
import java.beans.PropertyEditorSupport;
//Java2 imports
//import java.util.LinkedList;
//import java.util.Iterator;

import etomica.utility.LinkedList;
import etomica.utility.Iterator;

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
        LinkedList list = Simulation.instance.meterList();
        //count how many of all MeterAbstract objects are MeterFunction objects
        int nMeters = 1;

        for(Iterator iter=list.iterator(); iter.hasNext(); ) {
            if(iter.next() instanceof MeterFunction) nMeters++;
        }
        //loop again, save all the MeterFunction objects in the arrays
        meterNames = new String[nMeters];
        meterObjects = new Object[nMeters];
        meterNames[0] = "";
        meterObjects[0] = null;
        nMeters = 1;
        for(Iterator iter=list.iterator(); iter.hasNext(); ) {
            MeterAbstract meter = (MeterAbstract)iter.next();
            if(meter instanceof MeterFunction) {
                meterObjects[nMeters] = meter;
                meterNames[nMeters] = meterObjects[nMeters].toString();
                nMeters++;
            }
        }
    }
}