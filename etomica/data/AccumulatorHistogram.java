/*
 * History
 * Created on Aug 4, 2004 by kofke
 */
package etomica.data;

import etomica.Data;
import etomica.DataInfo;
import etomica.data.types.DataArithmetic;
import etomica.data.types.DataDoubleArray;
import etomica.units.Dimension;
import etomica.utility.Histogram;
import etomica.utility.HistogramSimple;

/**
 * Accumulator that keeps histogram of data.
 */
public class AccumulatorHistogram extends DataAccumulator {

    Histogram[] histogram = new Histogram[0];
    private Data data;
    private DataInfo binnedDataInfo;
    int nData, nDataMinus1;
    private Histogram.Factory histogramFactory;
    private int nBins;

    /**
     * Creates instance using HistogramSimple factory and specifying histograms
     * having 100 bins.
     */
    public AccumulatorHistogram(DataInfo info) {
        this(info, HistogramSimple.FACTORY);
    }

    public AccumulatorHistogram(DataInfo info, Histogram.Factory factory) {
        this(info, factory, 100);
    }

    public AccumulatorHistogram(DataInfo info, Histogram.Factory factory, int nBins) {
        this.nBins = nBins;
        this.binnedDataInfo = info;
        setNData(0);
        histogramFactory = factory;
    }

    /**
     * Adds each value in the given array to its own histogram. If the number of
     * values is different from that given in previous calls to the method, old
     * histogram data is discarded and new histograms are constructed (this
     * behavior can be modified by overriding the setNData method).
     */
    protected void addData(Data data) {
        DataArithmetic values = (DataArithmetic) data;
        if (values.getLength() != nData) {
            setNData(values.getLength());
        }
        for (int i = nDataMinus1; i >= 0; i--) {
            histogram[i].addValue(values.getValue(i));
        }
    }

    /**
     * Returns the set of histograms.
     */
    public Data getData() {
        if(nData == 1) {
            ((DataDoubleArray)((DataGroup)data).getData(0)).E(histogram[0].getHistogram());
            ((DataDoubleArray)((DataGroup)data).getData(1)).E(histogram[0].xValues());
            
        } else {
            for (int i = 0; i < nData; i++) {
                ((DataDoubleArray)((DataGroup)((DataGroup)data).getData(i)).getData(0)).E(histogram[i].getHistogram());
                ((DataDoubleArray)((DataGroup)((DataGroup)data).getData(i)).getData(1)).E(histogram[0].xValues());
            }
        }
        return data;
    }

    /**
     * Determines the number of histograms to be recorded, and is invoked when
     * the add method is called with an array argument of length different than
     * previously given. As implemented here, setNData(int) causes all
     * previously-recorded histograms to be discarded.
     * 
     * @param nData
     */
    protected void setNData(int nData) {
        this.nData = nData;
        nDataMinus1 = nData - 1;
        if(nData == 1) {
            DataDoubleArray dataH = new DataDoubleArray("Histogram", Dimension.NULL);
            DataDoubleArray dataBin = new DataDoubleArray(binnedDataInfo.getLabel(), binnedDataInfo.getDimension());
            data = new DataGroup("Histogram group", Dimension.NULL, new Data[] {dataH, dataBin});
        } else {
            DataGroup[] dataArray = new DataGroup[nData];
            for(int i=0; i<nData; i++) {
                DataDoubleArray dataH = new DataDoubleArray("Histogram", Dimension.NULL);
                DataDoubleArray dataBin = new DataDoubleArray(binnedDataInfo.getLabel(), binnedDataInfo.getDimension());
                dataArray[i] = new DataGroup("Histogram group", Dimension.NULL, new Data[] {dataH, dataBin});
            }
            data = new DataGroup("Group of Histogram groups", Dimension.NULL, dataArray);
        }
        histogram = new Histogram[nData];
        for (int i = 0; i < nData; i++)
            histogram[i] = histogramFactory.makeHistogram(nBins);
    }

    /**
     * @return Returns nBins, the number of bins in each histogram.
     */
    public int getNBins() {
        return nBins;
    }

    /**
     * @param bins
     *            Sets the number of bins in each histogram. Calls setNBins
     *            method of current histograms, which will discard data or
     *            modify themselves depending on how they are defined.
     */
    public void setNBins(int nBins) {
        this.nBins = nBins;
        for (int i = 0; i < nData; i++)
            histogram[i].setNBins(nBins);
    }

    public void reset() {
        for (int i = 0; i < nData; i++)
            histogram[i].reset();
    }

    /* (non-Javadoc)
     * @see etomica.DataSource#getDataInfo()
     */
    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }
}
