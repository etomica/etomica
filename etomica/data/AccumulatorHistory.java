/*
 * History
 * Created on Aug 11, 2004 by kofke
 */
package etomica.data;

import etomica.units.Undefined;
import etomica.data.types.CastToDoubleArray;
import etomica.data.types.DataArithmetic;
import etomica.data.types.DataArray;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataTable;
import etomica.units.Dimension;
import etomica.util.History;
import etomica.util.HistoryScrolling;

/**
 * Accumulator that keeps history of data.
 */
public class AccumulatorHistory extends DataAccumulator {

    History[] history = new History[0];
    DataSourceUniform xSource;
    private DataArray data;
    int nData;
    private History.Factory historyFactory;
    private int historyLength;

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
        xSource = new DataSourceUniform("Time", Undefined.DIMENSION, 1, historyLength, historyLength);
        this.historyLength = historyLength;
        historyFactory = factory;
    }

    /**
     * Returns caster that ensures accumulator will receive a DataDoubleArray.
     */
    public DataProcessor getDataCaster(DataInfo inputDataInfo) {
        if(inputDataInfo.getDataClass() == DataDoubleArray.class || 
                inputDataInfo.getDataClass() == DataTable.class) {
            return null;
        }
        return new CastToDoubleArray();
    }
    
    /**
     * Sets up data and histories, discarding any previous results.
     * 
     * @param nData
     */
    protected DataInfo processDataInfo(DataInfo inputDataInfo) {
        DataFactory factory = inputDataInfo.getDataFactory();
        if (factory instanceof DataDoubleArray.Factory) {
            nData = ((DataDoubleArray.Factory)inputDataInfo.getDataFactory()).getArrayLength();
        }
        else {
            nData = ((DataTable.Factory)inputDataInfo.getDataFactory()).getNRows()
                  * ((DataTable.Factory)inputDataInfo.getDataFactory()).getNColumns();
        }
        setupData(inputDataInfo);
        history = new History[nData];
        for (int i = 0; i < nData; i++) {
            history[i] = historyFactory.makeHistory(historyLength);
        }
        return data.getDataInfo();
    }
    
    /**
     * Adds each value in the given array to its own history. If the number of
     * values is different from that given in previous calls to the method, old
     * history data is discarded and new historys are constructed (this behavior
     * can be modified by overriding the setNData method).
     */
    protected void addData(Data data) {
        DataArithmetic values = (DataArithmetic)data;
        for (int i = nData-1; i >= 0; i--) {
            history[i].addValue(values.getValue(i));
        }
    }

    /**
     * Returns the set of histories.
     */
    public Data getData() {
        for (int i = 0; i < nData; i++) {
            DataFunction dataFunction = (DataFunction)data.getData(i);
            dataFunction.getYData().E(history[i].getHistory());
            dataFunction.getXData(0).E(history[i].getXSource().getData());
        }
        return data;
    }

    /**
     * Returns the number of histories times the length of each.
     */
    public int getDataLength() {
        return nData * historyLength;
    }

    public DataSource getXSource() {
        return history[0].getXSource();
    }

    /**
     * Constructs the Data objects used by this class.
     */
    private void setupData(DataInfo inputDataInfo) {
        DataDoubleArray dataBin = new DataDoubleArray("Time", Undefined.DIMENSION, historyLength);
        data = new DataArray(inputDataInfo.getLabel(), inputDataInfo.getDimension(), nData, DataFunction.getFactory(new DataDoubleArray[] {dataBin}));
    }
    
    /**
     * @return Returns historyLength, the number of bins in each history.
     */
    public int getHistoryLength() {
        return historyLength;
    }

    /**
     * @param bins
     *            Sets the number of bins in each history. Calls
     *            setHistoryLength method of current histories, which will
     *            discard data or modify themselves depending on how they are
     *            defined.
     */
    public void setHistoryLength(int historyLength) {
        this.historyLength = historyLength;
        xSource.setNValues(historyLength);
        xSource.setXMax(historyLength);
        for (int i = 0; i < nData; i++)
            history[i].setHistoryLength(historyLength);
    }

    public void reset() {
        for (int i = 0; i < nData; i++)
            history[i].reset();
    }
    
    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }
}
