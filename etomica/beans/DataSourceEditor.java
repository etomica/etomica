package etomica;
import java.beans.PropertyEditorSupport;

/**
 * Editor for selection of data source objects.
 * Provides a combo-box selector of all objects added to the simulation 
 * that implement the DataSource interface.
 */
public class DataSourceEditor extends PropertyEditorSupport implements java.io.Serializable {
    
    public static String version() {return "01.04.28";}
    
    Object[] dataSourceObjects;
    String[] dataSourceNames;
    
    public DataSourceEditor() {
        super();
        makeDataSourceLists();
    }

    public String getAsText() {
        makeDataSourceLists();
        if(dataSourceNames.length == 0 || getValue() == null) return "null";
        return getValue().toString();
    }
    
    public void setAsText(String s) {
        makeDataSourceLists();
        for(int i=0; i<dataSourceNames.length; i++) {
            if(dataSourceNames[i].equals(s)) {
                setValue(dataSourceObjects[i]);
                firePropertyChange();
                return;
            }
        }
    }

    public String[] getTags() {return dataSourceNames;}

    protected void makeDataSourceLists() {
        java.util.LinkedList dataSourceList = etomica.Simulation.instance.allElements();
        //count number of DataSource objects
        int count=0;
        for(java.util.Iterator iter=dataSourceList.iterator(); iter.hasNext(); ) {
            Object obj = iter.next();
            if(obj instanceof DataSource) count++;
        }
        //collect them together
        dataSourceNames = new String[count];
        dataSourceObjects = new Object[count];
        int i=0;
        for(java.util.Iterator iter=dataSourceList.iterator(); iter.hasNext(); ) {
            Object obj = iter.next();
            if(obj instanceof DataSource) {
                dataSourceObjects[i] = (DataSource)obj;
                dataSourceNames[i] = dataSourceObjects[i].toString();
                i++;
            }
        }
    }
}//end of class