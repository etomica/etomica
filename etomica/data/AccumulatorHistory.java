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
import etomica.data.types.DataFunction;
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
            ((DataFunction)data).E(history[0].getHistory());
            //TODO put this out of its misery
            ((DataFunction)data).getTData().E(((DataDoubleArray)history[0].getXSource().getData()).toArray());
            
        } else {
            for (int i = 0; i < nData; i++) {
                ((DataFunction)((DataGroup)data).getData(i)).E(history[i].getHistory());
                //TODO this too
                ((DataFunction)((DataGroup)data).getData(i)).getTData().E(((DataDoubleArray)history[i].getXSource().getData()).toArray());
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
            data = new DataDoubleArray(new DataInfo("History", Dimension.NULL));
        } else {
            DataInfo dataInfo = new DataInfo("History", Dimension.NULL);
            DataDoubleArray[] dataArray = new DataDoubleArray[nData];
            for (int i = 0; i < nData; i++) {
                dataArray[i] = new DataDoubleArray(dataInfo);
            }
            data = new DataGroup(new DataInfo("History Group", Dimension.NULL),
                    dataArray);
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
