/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.units.dimensions.Null;
import etomica.util.Arrays;


/**
 * Collects one or more data sources into a single DataSource. Data returned
 * by the DataSourceGroup is a DataGroup formed from the Data given
 * by the component DataSource instances.
 *
 * @author David Kofke
 *
 */
public class DataSourceGroup implements IEtomicaDataSource, java.io.Serializable {

    /**
     * Forms a DataSourceGroup that contains no data sources, and which will give
     * an empty DataGroup for getData.  Must populate group with addDataSource method. 
     */
    public DataSourceGroup() {
        data = new DataGroup(new IData[0]);
        dataInfo = new DataInfoGroup("Data Group", Null.DIMENSION, new IEtomicaDataInfo[0]);
        dataSources = new IEtomicaDataSource[0];
        tag = new DataTag();
        dataInfo.addTag(tag);
    }
    
    /**
     * Forms a DataGroup containing the given data sources.  Given
     * array is copied to another array internally.
     */
    public DataSourceGroup(IEtomicaDataSource[] sources) {
        this();
        if(sources != null) {
            for(int i=0; i<sources.length; i++) {
                addDataSource(sources[i]);
            }
        }
    }
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }

    /**
     * Returns a DataGroup that is formed from the Data objects returned by
     * the data sources in this instance.  Dimension of returned data is that of
     * the Data it holds, if they are all the same; otherwise it is Dimension.UNDEFINED.
     */
    public IData getData() {
        boolean rebuildData = false;
        //generate data from sources, check that all returned instances are the same as before
        //if a new data instance is given, make a new DataGroup that uses it instead of the previous data
        for(int i=0; i<dataSources.length; i++) {
            if(latestData[i] != data.getData(i)) {
                rebuildData = true;
            }
        }
        if(rebuildData) {
            dataInfo = new DataInfoGroup("Data Group", Null.DIMENSION, getSubDataInfo());
            dataInfo.addTag(tag);
            data = new DataGroup(latestData);
        }
        return data;
    }
    
    /**
     * Adds the given DataSource to those held by the group.
     */
    public void addDataSource(IEtomicaDataSource newSource) {
        dataSources = (IEtomicaDataSource[])Arrays.addObject(dataSources, newSource);
        latestData = (IData[])Arrays.addObject(latestData, newSource.getData());
        dataInfo = new DataInfoGroup("Data Group", Null.DIMENSION, getSubDataInfo());
        data = new DataGroup(latestData);
    }
    
    protected IEtomicaDataInfo[] getSubDataInfo() {
        IEtomicaDataInfo[] subDataInfo = new IEtomicaDataInfo[dataSources.length];
        for (int i=0; i<subDataInfo.length; i++) {
            subDataInfo[i] = dataSources[i].getDataInfo();
        }
        return subDataInfo;
    }
    
    private static final long serialVersionUID = 1L;
    private IEtomicaDataSource[] dataSources = new IEtomicaDataSource[0];
    private IData[] latestData = new IData[0];
    private DataGroup data;
    private IEtomicaDataInfo dataInfo;
    protected final DataTag tag;
}
