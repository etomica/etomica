package etomica.data;

import etomica.data.types.DataGroup;
import etomica.utility.Arrays;


/**
 * Collects one or more data sources into a single DataSource. Data returned
 * by the DataSourceGroup is a DataGroup formed from the Data given
 * by the component DataSource instances.
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Apr 8, 2005 by kofke
 */
public class DataSourceGroup implements DataSource, java.io.Serializable {

    /**
     * Forms a DataSourceGroup that contains no data sources, and which will give
     * an empty DataGroup for getData.  Must populate group with addDataSource method. 
     */
    public DataSourceGroup() {
        data = new DataGroup("Data Group", new Data[0]);
        dataSources = new DataSource[0];
    }
    
    /**
     * Forms a DataGroup containing the given data sources.  Given
     * array is copied to another array internally.
     */
    public DataSourceGroup(DataSource[] sources) {
        this();
        if(sources != null) {
            for(int i=0; i<sources.length; i++) {
                addDataSource(sources[i]);
            }
        }
    }
    
    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }

    /**
     * Returns a DataGroup that is formed from the Data objects returned by
     * the data sources in this instance.  Dimension of returned data is that of
     * the Data it holds, if they are all the same; otherwise it is Dimension.UNDEFINED.
     */
    public Data getData() {
        boolean rebuildData = false;
        //generate data from sources, check that all returned instances are the same as before
        //if a new data instance is given, make a new DataGroup that uses it instead of the previous data
        for(int i=0; i<dataSources.length; i++) {
            latestData[i] = dataSources[i].getData();
            if(latestData[i] != data.getData(i)) {
                rebuildData = true;
            }
        }
        if(rebuildData) {
            data = new DataGroup(data.getDataInfo().getLabel(), latestData);
        }
    	    return data;
    }
    
    /**
     * Adds the given DataSource to those held by the group.
     */
    public void addDataSource(DataSource newSource) {
        dataSources = (DataSource[])Arrays.addObject(dataSources, newSource);
        latestData = (Data[])Arrays.addObject(latestData, newSource.getData());
        data = new DataGroup("Data Group", latestData);
    }

    private DataSource[] dataSources = new DataSource[0];
    private Data[] latestData = new Data[0];
    private DataGroup data;
}
