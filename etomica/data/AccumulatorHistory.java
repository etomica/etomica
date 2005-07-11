/*
 * History
 * Created on Aug 11, 2004 by kofke
 */
package etomica.data;

import etomica.Data;
import etomica.DataInfo;
import etomica.DataSource;
import etomica.data.types.DataArithmetic;
import etomica.data.types.DataDoubleArray;
import etomica.units.Dimension;
import etomica.utility.History;
import etomica.utility.HistoryScrolling;

/**
 * Accumulator that keeps history of data.
 */
public class AccumulatorHistory extends DataAccumulator {

    History[] history = new History[0];
    DataSourceUniform xSource;
    private Data data;
    int nData, nDataMinus1;
    private DataInfo binnedDataInfo;
    private History.Factory historyFactory;
    private int historyLength;

    /**
     * Creates instance using HistorySimple factory and specifying historys
     * having 100 bins.
     */
    public AccumulatorHistory(DataInfo info) {
        this(info, HistoryScrolling.FACTORY);
    }

    public AccumulatorHistory(DataInfo info, History.Factory factory) {
        this(info, factory, 100);
    }

    public AccumulatorHistory(DataInfo info, History.Factory factory, int historyLength) {
        super();
        binnedDataInfo = info;
        xSource = new DataSourceUniform(Dimension.UNDEFINED, 1, historyLength, historyLength);
        this.historyLength = historyLength;
        historyFactory = factory;
        setNData(1);
    }

    /**
     * Adds each value in the given array to its own history. If the number of
     * values is different from that given in previous calls to the method, old
     * history data is discarded and new historys are constructed (this behavior
     * can be modified by overriding the setNData method).
     */
    protected void addData(Data data) {
        DataArithmetic values = (DataArithmetic)data;
        if (values.getLength() != nData) {
            setNData(values.getLength());
        }
        for (int i = nDataMinus1; i >= 0; i--) {
            history[i].addValue(values.getValue(i));
        }
    }

    /**
     * Returns the set of histograms.
     */
    public Data getData() {
        if(nData == 1) {
            ((DataDoubleArray)((DataGroup)data).getData(0)).E(history[0].getHistory());
            ((DataDoubleArray)((DataGroup)data).getData(1)).E(history[0].getXSource().getData());
           
        } else {
            for (int i = 0; i < nData; i++) {
                ((DataDoubleArray)((DataGroup)((DataGroup)data).getData(i)).getData(0)).E(history[i].getHistory());
                ((DataDoubleArray)((DataGroup)((DataGroup)data).getData(i)).getData(1)).E(history[0].getXSource().getData());
            }
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
     * Determines the number of historys to be recorded, and is invoked when the
     * add method is called with an array argument of length different than
     * previously given. As implemented here, setNData(int) causes all
     * previously-recorded histories to be discarded.
     * 
     * @param nData
     */
    protected void setNData(int nData) {
        this.nData = nData;
        nDataMinus1 = nData - 1;

        if (nData == 1) {
            DataDoubleArray dataH = new DataDoubleArray(binnedDataInfo.getLabel()+" History", Dimension.NULL);
            DataDoubleArray dataBin = new DataDoubleArray("Time", Dimension.UNDEFINED);
            dataH.setLength(0);
            dataBin.setLength(0);
            data = new DataGroup("History group", Dimension.NULL, new Data[] {dataH, dataBin});

        } else {
            DataGroup[] dataArray = new DataGroup[nData];
            for(int i=0; i<nData; i++) {
                DataDoubleArray dataH = new DataDoubleArray(binnedDataInfo.getLabel()+" History", Dimension.NULL);
                DataDoubleArray dataBin = new DataDoubleArray("Time", Dimension.UNDEFINED);
                dataH.setLength(0);
                dataBin.setLength(0);
                dataArray[i] = new DataGroup("History group", Dimension.NULL, new Data[] {dataH, dataBin});
            }
            data = new DataGroup("Group of History groups", Dimension.NULL, dataArray);
        }

        history = new History[nData];
        for (int i = 0; i < nData; i++) {
            history[i] = historyFactory.makeHistory(historyLength);
        }
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
