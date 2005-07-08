package etomica.data;

import etomica.Data;
import etomica.DataInfo;
import etomica.DataSource;
import etomica.units.Dimension;
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
        data = new DataGroup(new DataInfo("Data Group", Dimension.NULL), new Data[0]);
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
     * the data sources in this instance.
     */
    public Data getData() {
        for(int i=0; i<dataSources.length; i++) {
            data.data[i] = dataSources[i].getData();
        }
    	    return data;
    }
    
    /**
     * Adds the given DataSource to those held by the group.  The new data source
     * must give data having the same physical dimensions (as indicated by its
     * getDimension()) as those already added (if any).
     */
    public void addDataSource(DataSource newSource) {
        Dimension newDimension = newSource.getDataInfo().getDimension();
        if(data.getNData() > 0 && (newDimension != data.getDataInfo().getDimension())) {
            newDimension = Dimension.UNDEFINED; //TODO maybe make a Dimension.MIXED 
        }
        dataSources = (DataSource[])Arrays.addObject(dataSources, newSource);
        data = new DataGroup(new DataInfo("Data Group", newDimension), new Data[dataSources.length]);
    }

    private DataSource[] dataSources;
    private DataGroup data;
}
