/*
 * History
 * Created on Aug 11, 2004 by kofke
 */
package etomica.data;

import etomica.data.types.CastToDoubleArray;
import etomica.data.types.DataArithmetic;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataGroup;
import etomica.data.types.DataTable;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.data.types.DataTable.DataInfoTable;
import etomica.units.Quantity;
import etomica.util.History;
import etomica.util.HistoryScrolling;

/**
 * Accumulator that keeps history of data.
 */
public class AccumulatorHistory extends DataAccumulator {

    /**
     * Creates instance using HistorySimple factory and specifying historys
     * having 100 bins.
     */
    public AccumulatorHistory() {
        this(HistoryScrolling.FACTORY);
    }

    public AccumulatorHistory(History.Factory factory) {
        this(factory, 100);
    }

    public AccumulatorHistory(History.Factory factory, int historyLength) {
        super();
        this.historyLength = historyLength;
        historyFactory = factory;
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
     * Returns caster that ensures accumulator will receive a DataDoubleArray.
     */
    public DataProcessor getDataCaster(DataInfo newInputDataInfo) {
        if(newInputDataInfo instanceof DataInfoDoubleArray || 
                newInputDataInfo instanceof DataInfoTable) {
            return null;
        }
        return new CastToDoubleArray();
    }
    
    /**
     * Sets up data and histories, discarding any previous results.
     * 
     * @param nData
     */
    protected DataInfo processDataInfo(DataInfo newInputDataInfo) {
        inputDataInfo = newInputDataInfo;
        if (inputDataInfo instanceof DataInfoDoubleArray) {
            nData = ((DataInfoDoubleArray)inputDataInfo).getArrayLength();
        }
        else {
            //if it's not a DataDoubleArray, it must be a DataTable
            nData = ((DataInfoTable)inputDataInfo).getNRows()
                  * ((DataInfoTable)inputDataInfo).getNDataInfo();
        }
        history = new History[nData];
        for (int i = 0; i < nData; i++) {
            history[i] = historyFactory.makeHistory(historyLength);
        }
        setupData();
        return dataInfo;
    }
    
    /**
     * Adds each value in the given array to its own history. If the number of
     * values is different from that given in previous calls to the method, old
     * history data is discarded and new historys are constructed (this behavior
     * can be modified by overriding the setNData method).
     */
    protected void addData(Data newData) {
        DataArithmetic values = (DataArithmetic)newData;
        for (int i = nData-1; i >= 0; i--) {
            history[i].addValue(timeDataSource.getDataAsScalar(), values.getValue(i));
        }
    }

    /**
     * Returns the set of histories.
     */
    public Data getData() {
        // check to see if current data functions are the right length
        // histogram might change the number of bins on its own.
        boolean success = true;
        for (int i=0; i<nData; i++) {
            DataFunction dataFunction = (DataFunction)data.getData(i);
            // covertly call getHistogram() here.  We have to call it anyway
            // so that the array data gets updated.
            if (dataFunction.getData() != history[i].getHistory() ||
                    dataFunction.getXData(0).getData() != history[i].getXValues()) {
                success = false;
            }
        }

        if (!success) {
            DataFunction[] dataFunctions = new DataFunction[nData];
            DataInfoFunction[] dataInfoFunctions = new DataInfoFunction[nData];
            // attempt to re-use old DataFunctions
            for (int i=0; i<nData; i++) {
                DataFunction dataFunction = (DataFunction)data.getData(i);
                DataInfoFunction dataInfoFunction = (DataInfoFunction)((DataInfoGroup)dataInfo).getSubDataInfo(i);
                if (dataFunction.getData() == history[i].getHistory() &&
                        dataFunction.getXData(0).getData() == history[i].getXValues()) {
                    dataFunctions[i] = dataFunction;
                    dataInfoFunctions[i] = dataInfoFunction;
                }
                else {
                    int iHistoryLength = history[i].getHistoryLength();
                    dataInfoFunctions[i] = new DataInfoFunction(inputDataInfo.getLabel(), inputDataInfo.getDimension(), 
                            new DataInfoDoubleArray[]{new DataInfoDoubleArray(timeDataSource.getDataInfo().getLabel(),
                                    timeDataSource.getDataInfo().getDimension(), new int[]{iHistoryLength})});
                    DataDoubleArray xData = new DataDoubleArray(new int[]{history[i].getHistoryLength()},history[i].getXValues());
                    dataFunctions[i] = new DataFunction(new DataDoubleArray[]{xData}, history[i].getHistory());
                }
            }
            // creating a new data instance might confuse downstream data sinks, but
            // we have little choice and they should deal.
            data = new DataGroup(dataFunctions);
            dataInfo = new DataInfoGroup(inputDataInfo.getLabel(), inputDataInfo.getDimension(), dataInfoFunctions);
        }
            
        return data;
    }

    /**
     * Returns the number of histories times the length of each.
     */
    public int getDataLength() {
        return nData * historyLength;
    }

    /**
     * Constructs the Data objects used by this class.
     */
    private void setupData() {
        DataFunction[] dataFunctions = new DataFunction[nData];
        DataInfoFunction[] dataInfoFunctions = new DataInfoFunction[nData];
        for (int i=0; i<nData; i++) {
            dataInfoFunctions[i] = new DataInfoFunction(inputDataInfo.getLabel(), inputDataInfo.getDimension(), 
                    new DataInfoDoubleArray[]{new DataInfoDoubleArray(timeDataSource.getDataInfo().getLabel(),
                            timeDataSource.getDataInfo().getDimension(), new int[]{history[i].getHistoryLength()})});
            DataDoubleArray xData = new DataDoubleArray(new int[]{history[i].getHistoryLength()},history[i].getXValues());
            dataFunctions[i] = new DataFunction(new DataDoubleArray[]{xData}, history[i].getHistory());
        }
        data = new DataGroup(dataFunctions);
        dataInfo = new DataInfoGroup(inputDataInfo.getLabel(), inputDataInfo.getDimension(), dataInfoFunctions);
    }
    
    /**
     * @return Returns historyLength, the number of bins in each history.
     */
    public int getHistoryLength() {
        return historyLength;
    }

    /**
     * @param historyLength
     *            Sets the number of bins in each history. Calls
     *            setHistoryLength method of current histories, which will
     *            discard data or modify themselves depending on how they are
     *            defined.
     */
    public void setHistoryLength(int historyLength) {
        this.historyLength = historyLength;
        for (int i = 0; i < nData; i++)
            history[i].setHistoryLength(historyLength);
    }

    public void reset() {
        for (int i = 0; i < nData; i++)
            history[i].reset();
    }
    
    public DataInfo getDataInfo() {
        return dataInfo;
    }
    
    protected History[] history = new History[0];
    private DataGroup data;
    int nData;
    private History.Factory historyFactory;
    private int historyLength;
    private DataSourceScalar timeDataSource;
    private DataInfo inputDataInfo;

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
