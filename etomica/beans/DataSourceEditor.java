package etomica;
import java.beans.PropertyEditorSupport;
//Java2 imports
//import java.util.LinkedList;
//import java.util.Iterator;

import etomica.utility.LinkedList;
import etomica.utility.Iterator;

/**
 * Editor for selection of data source objects.
 * Provides a combo-box selector of all objects added to the simulation 
 * that implement the DataSource interface.
 */
public class DataSourceEditor extends PropertyEditorSupport implements java.io.Serializable {
    
    public static String version() {return "01.05.14";}
    
    Object[] dataSourceObjects;
    String[] dataSourceNames;
    String valueAsText;
    
    public DataSourceEditor() {
        super();
        makeDataSourceList();
    }

    public String getAsText() {
//        return valueAsText;
        makeDataSourceList();
        if(dataSourceNames.length == 0 || getValue() == null) return "null";
        return getValue().toString();
    }
    
    public void setAsText(String s) {
        makeDataSourceList();
        for(int i=0; i<dataSourceNames.length; i++) {
            if(dataSourceNames[i].equals(s)) {
                int idx = s.indexOf(':');
                if(idx == -1) {  //not a wrapper
                    setValue(dataSourceObjects[i]);
                }
                else { //wrapper
                    String sourceText = s.substring(idx+1);
                    DataSource.Wrapper wrapper = (DataSource.Wrapper)dataSourceObjects[i];
                    setValue(wrapper.getDataSource(sourceText));
                }
                firePropertyChange();
                valueAsText = s;
                return;
            }
        }
    }

    public String[] getTags() {return dataSourceNames;}

    protected void makeDataSourceList() {
        LinkedList dataSourceList = etomica.Simulation.instance.allElements();
        //count number of DataSource objects
        int count=1;
        for(Iterator iter=dataSourceList.iterator(); iter.hasNext(); ) {
            Object obj = iter.next();
            if(obj instanceof DataSource) count++;
            if(obj instanceof DataSource.Wrapper) 
                count += ((DataSource.Wrapper)obj).getSourcesAsText().length;
        }
//        if(count > 2) count++; //more than one source; will add "More..." entry
        //collect them together
        dataSourceNames = new String[count];
        dataSourceObjects = new Object[count];
        dataSourceNames[0] = "null";
        dataSourceObjects[0] = null;
        int i=1;
        for(Iterator iter=dataSourceList.iterator(); iter.hasNext(); ) {
            Object obj = iter.next();
            if(obj instanceof DataSource) {
                dataSourceObjects[i] = (DataSource)obj;
                dataSourceNames[i] = dataSourceObjects[i].toString();
                i++;
            }
            if(obj instanceof DataSource.Wrapper) {
                DataSource.Wrapper wrapper = (DataSource.Wrapper)obj;
                String prefix = wrapper.toString() + ':';
                String[] names = wrapper.getSourcesAsText();
                for(int j=0; j<names.length; j++) {
                    dataSourceObjects[i] = wrapper;
                    dataSourceNames[i] = prefix+names[j];
                    i++;
                }
            }//end if
        }//end of for
//        if(count > 2) {
//            dataSourceNames[i] = "More...";
//        }
    }//end of makeDataSourceList
}//end of class