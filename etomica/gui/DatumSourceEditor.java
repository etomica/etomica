package etomica;
import java.beans.PropertyEditorSupport;
//Java2 imports
//import java.util.LinkedList;
//import java.util.Iterator;

import etomica.utility.LinkedList;
import etomica.utility.Iterator;

/**
 * Editor for selection of datum source objects.
 * Provides a combo-box selector of all objects added to the simulation 
 * that implement the DatumSource interface.
 */
public class DatumSourceEditor extends PropertyEditorSupport implements java.io.Serializable {
    
    public static String version() {return "01.05.14";}
    
    Object[] datumSourceObjects;
    String[] datumSourceNames;
    String valueAsText;
    
    public DatumSourceEditor() {
        super();
        makeDatumSourceList();
    }

    public String getAsText() {
//        return valueAsText;
        makeDatumSourceList();
        if(datumSourceNames.length == 0 || getValue() == null) return "null";
        return getValue().toString();
    }
    
    public void setAsText(String s) {
        makeDatumSourceList();
        for(int i=0; i<datumSourceNames.length; i++) {
            if(datumSourceNames[i].equals(s)) {
                int idx = s.indexOf(':');
//                if(idx == -1) {  //not a wrapper
                    setValue(datumSourceObjects[i]);
/*                }
                else { //wrapper
                    String sourceText = s.substring(idx+1);
                    DatumSource.Wrapper wrapper = (DatumSource.Wrapper)datumSourceObjects[i];
                    setValue(wrapper.getDatumSource(sourceText));
                }*/
                firePropertyChange();
                valueAsText = s;
                return;
            }
        }
    }

    public String[] getTags() {return datumSourceNames;}

    protected void makeDatumSourceList() {
        LinkedList datumSourceList = etomica.Simulation.instance.allElements();
        //count number of DatumSource objects
        int count=1;
        for(Iterator iter=datumSourceList.iterator(); iter.hasNext(); ) {
            Object obj = iter.next();
            if(obj instanceof DatumSource) count++;
//            if(obj instanceof DatumSource.Wrapper) 
//                count += ((DatumSource.Wrapper)obj).getSourcesAsText().length;
        }
//        if(count > 2) count++; //more than one source; will add "More..." entry
        //collect them together
        datumSourceNames = new String[count];
        datumSourceObjects = new Object[count];
        datumSourceNames[0] = "null";
        datumSourceObjects[0] = null;
        int i=1;
        for(Iterator iter=datumSourceList.iterator(); iter.hasNext(); ) {
            Object obj = iter.next();
            if(obj instanceof DatumSource) {
                datumSourceObjects[i] = (DatumSource)obj;
                datumSourceNames[i] = datumSourceObjects[i].toString();
                i++;
            }
/*            if(obj instanceof DatumSource.Wrapper) {
                DatumSource.Wrapper wrapper = (DatumSource.Wrapper)obj;
                String prefix = wrapper.toString() + ':';
                String[] names = wrapper.getSourcesAsText();
                for(int j=0; j<names.length; j++) {
                    datumSourceObjects[i] = wrapper;
                    datumSourceNames[i] = prefix+names[j];
                    i++;
                }
            }//end if  */
        }//end of for
//        if(count > 2) {
//            datumSourceNames[i] = "More...";
//        }
    }//end of makeDatumSourceList
}//end of class