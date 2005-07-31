/*
 * History
 * Created on Aug 4, 2004 by kofke
 */
package etomica.data;

import etomica.Data;
import etomica.DataInfo;
import etomica.data.types.CastToDoubleArray;
import etomica.data.types.DataArithmetic;
import etomica.data.types.DataArray;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.units.Dimension;
import etomica.utility.Histogram;
import etomica.utility.HistogramSimple;

/**
 * Accumulator that keeps histogram of data.
 */
public class AccumulatorHistogram extends DataAccumulator {

    /**
     * Creates instance using HistogramSimple factory and specifying histograms
     * having 100 bins.
     */
    public AccumulatorHistogram() {
        this(HistogramSimple.FACTORY);
    }

    public AccumulatorHistogram(Histogram.Factory factory) {
        this(factory, 100);
    }

    public AccumulatorHistogram(Histogram.Factory factory, int nBins) {
        this.nBins = nBins;
        histogramFactory = factory;
    }

    /**
     * Adds each value in the given array to its own histogram. If the number of
     * values is different from that given in previous calls to the method, old
     * histogram data is discarded and new histograms are constructed (this
     * behavior can be modified by overriding the setNData method).
     */
    protected void addData(Data inputData) {
        DataArithmetic values = (DataArithmetic)inputData;
        for (int i = nData-1; i >= 0; i--) {
            histogram[i].addValue(values.getValue(i));
        }
    }

    /**
     * Returns the set of histograms.
     */
    public Data getData() {
        for (int i = 0; i < nData; i++) {
            DataFunction dataFunction = (DataFunction)data.getData(i);
            dataFunction.getYData().E(histogram[i].getHistogram());
            dataFunction.getXData(0).E(histogram[i].xValues());
        }
        return data;
    }
    
    /* (non-Javadoc)
     * @see etomica.data.DataProcessor#makeOutputDataInfo(etomica.DataInfo)
     */
    /**
     * Sets up data and histograms, discarding any previous results.
     * 
     * @param nData
     */
    protected DataInfo processDataInfo(DataInfo inputDataInfo) {
        binnedDataInfo = inputDataInfo;
        nData = ((DataDoubleArray.Factory)inputDataInfo.getDataFactory()).getArrayLength();
        setupData();
        histogram = new Histogram[nData];
        for (int i = 0; i < nData; i++) {
            histogram[i] = histogramFactory.makeHistogram(nBins);
        }
        return data.getDataInfo();
    }
    
    /**
     * Returns caster that ensures accumulator will receive a DataDoubleArray.
     */
    public DataProcessor getDataCaster(DataInfo inputDataInfo) {
        if(inputDataInfo.getClass() == DataDoubleArray.class) {
            return null;
        }
        return new CastToDoubleArray();
    }

    /**
     * Constructs the Data objects used by this class.
     */
    private void setupData() {
        DataDoubleArray dataBin = new DataDoubleArray(binnedDataInfo.getLabel(), binnedDataInfo.getDimension(), nBins);
//        dataH = new DataDoubleArray("Histogram", Dimension.NULL, new int[] {nBins, nData});
//        data = new DataFunction(new DataDoubleArray[] {dataBin}, dataH);
        data = new DataArray("Histogram", Dimension.NULL, nData, DataFunction.getFactory(new DataDoubleArray[] {dataBin}));
    }
    
    /**
     * @return Returns nBins, the number of bins in each histogram.
     */
    public int getNBins() {
        return nBins;
    }

    /**
     * Sets the number of bins in each histogram. Calls setNBins method of
     * current histograms, which will discard data or modify themselves
     * depending on how they are defined.
     */
    public void setNBins(int nBins) {
        this.nBins = nBins;
        for (int i = 0; i < nData; i++) {
            histogram[i].setNBins(nBins);
        }
        setupData();
        if(dataSink != null) {
            dataSink.putDataInfo(data.getDataInfo());
        }
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
    
    Histogram[] histogram = new Histogram[0];
    private DataArray data;
    private DataInfo binnedDataInfo;
    int nData;
    private Histogram.Factory histogramFactory;
    private int nBins;
//    private DataDoubleArray dataH, dataBin;


}
