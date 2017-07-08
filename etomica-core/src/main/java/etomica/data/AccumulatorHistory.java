/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.units.dimensions.Quantity;
import etomica.data.history.History;
import etomica.data.history.HistoryScrolling;

/**
 * Accumulator that keeps history of data.
 */
public class AccumulatorHistory extends DataAccumulator {

    /**
     * Creates instance using HistorySimple factory and specifying historys
     * having 100 bins.
     */
    public AccumulatorHistory() {
        this(new HistoryScrolling());
    }

    public AccumulatorHistory(History history) {
        this.history = history;
        setTimeDataSource(new DataSourceCount());
    }

    /**
     * Sets the DataSource used for the "time" component of the history.  Each
     * time addData is called, the time DataSource's getData will be called and
     * the returned scalar will be taken as the current time.  By default, the 
     * time variable is taken to be the number of times data was added to the 
     * history.
     */
    public void setTimeDataSource(DataSourceScalar newTimeDataSource) {
        timeDataSource = newTimeDataSource;
        if (dataInfo != null) {
            putDataInfo(inputDataInfo);
        }
    }
    
    /**
     * Returns the DataSource used for the "time" component of the history.
     */
    public DataSourceScalar getTimeDataSource() {
        return timeDataSource;
    }
    
    /**
     * Returns null.  AccumulatorHistogram can take an type of Data.
     */
    public DataPipe getDataCaster(IEtomicaDataInfo newInputDataInfo) {
        return null;
    }
    
    /**
     * Sets up data and histories, discarding any previous results.
     * 
     * @param newInputDataInfo
     */
    protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo newInputDataInfo) {
        if (newInputDataInfo.getLength() != 1) {
            throw new IllegalArgumentException("AccumulatorHistory only handles single-value data");
        }
        inputDataInfo = newInputDataInfo;
        setupData();
        return dataInfo;
    }
    
    /**
     * Adds each value in the given array to its own history. If the number of
     * values is different from that given in previous calls to the method, old
     * history data is discarded and new historys are constructed (this behavior
     * can be modified by overriding the setNData method).
     */
    protected boolean addData(IData newData) {
        return history.addValue(timeDataSource.getDataAsScalar(), newData.getValue(0));
    }

    /**
     * Returns the set of histories.
     */
    public IData getData() {
        // check to see if current data functions are the right length
        // histogram might change the number of bins on its own.
        // covertly call getHistogram() here.  We have to call it anyway
        // so that the array data gets updated.
        if (data.getData() != history.getHistory() ||
            xDataSources.getIndependentData(0).getData() != history.getXValues()) {
            setupData();
        }
        
        return data;
    }

    /**
     * Constructs the Data objects used by this class.
     */
    private void setupData() {
        data = new DataFunction(new int[]{history.getHistoryLength()}, history.getHistory());
        xDataSources = new DataSourceIndependentSimple(history.getXValues(), 
                    new DataInfoDoubleArray(timeDataSource.getDataInfo().getLabel(),
                    timeDataSource.getDataInfo().getDimension(), new int[]{history.getHistoryLength()}));
        dataInfo = new DataInfoFunction(inputDataInfo.getLabel(), 
                    inputDataInfo.getDimension(), xDataSources);
        dataInfo.addTags(inputDataInfo.getTags());
        dataInfo.addTag(tag);
        if (dataSink != null) {
            dataSink.putDataInfo(dataInfo);
        }
    }
    
    public History getHistory() {
        return history;
    }
    
    public void setHistory(History newHistory) {
        history = newHistory;
    }
    
    public void reset() {
        history.reset();
    }
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    protected History history;
    protected DataSourceIndependentSimple xDataSources;
    private DataFunction data;
    protected int nData;
    private DataSourceScalar timeDataSource;
    private IEtomicaDataInfo inputDataInfo;

    /**
     * Simple DataSource to use as a default time DataSource.  It just returns
     * the number of times it's been called.
     */
    protected static class DataSourceCount extends DataSourceScalar {

        public DataSourceCount() {
            super("Count",Quantity.DIMENSION);
        }
        
        public double getDataAsScalar() {
            return count++;
        }

        private int count = 0;
    }
}
