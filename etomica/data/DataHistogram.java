/*
 * History
 * Created on Aug 4, 2004 by kofke
 */
package etomica.data;

import etomica.data.types.CastToDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.units.Null;
import etomica.util.DoubleRange;
import etomica.util.Histogram;
import etomica.util.HistogramSimple;

/**
 * DataProcessor that creates a histogram from each piece of data that comes in.
 * A single Histogram object is kept, but is reset for each new data.
 * <p>
 * Input Data must implement DataArithmetic.
 */
public class DataHistogram extends DataProcessor {

    /**
     * Creates instance using HistogramSimple factory and specifying histograms
     * having the given number of bins and the given range.
     */
    public DataHistogram(int nBins, DoubleRange range) {
        this(new HistogramSimple.Factory(nBins, range), nBins);
    }

    /**
     * Creates instance using given histogram factory with default nBins of 100.
     */
    public DataHistogram(Histogram.Factory factory) {
        this(factory, 100);
    }

    /**
     * Creates instance using the given histogram factory making histograms having
     * the given number of bins.
     */
    public DataHistogram(Histogram.Factory factory, int nBins) {
        this.nBins = nBins;
        histogramFactory = factory;
    }

    /**
     * Adds each value in the given Data to a single histogram.
     */
    protected Data processData(Data inputData) {
        histogram.reset();
        for (int i = 0; i <nData; i++) {
            histogram.addValue(inputData.getValue(i));
        }

        // check to see if current data function is the right length
        // histogram might change the number of bins on its own.
        boolean success = true;
        // covertly call getHistogram() here.  We have to call it anyway
        // so that the array data gets updated.
        if (data.getData() != histogram.getHistogram() ||
                xDataSource.getIndependentData(0).getData() != histogram.xValues()) {
            success = false;
        }
        
        if (!success) {
            DataInfoDoubleArray independentInfo = new DataInfoDoubleArray(binnedDataInfo.getLabel(),binnedDataInfo.getDimension(),new int[]{histogram.getNBins()});
            data = new DataFunction(new int[]{histogram.getNBins()}, histogram.getHistogram());
            xDataSource = new DataSourceIndependentSimple(histogram.xValues(), independentInfo);
            dataInfo = new DataInfoFunction(binnedDataInfo.getLabel()+" Histogram", Null.DIMENSION, xDataSource);

            // if we have have our own DataSink, we need to notify it that something changed.
            if(dataSink != null) {
                dataSink.putDataInfo(dataInfo);
            }
        }
        
        return data;
    }
    
    /**
     * Sets up the histogram, discarding any previous results.
     */
    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        binnedDataInfo = inputDataInfo;
        nData = ((DataInfoDoubleArray)inputDataInfo).getArrayLength();
        histogram = histogramFactory.makeHistogram();
        histogram.setNBins(nBins);
        setupData();
        return dataInfo;
    }
    
    /**
     * Returns caster that ensures accumulator will receive a DataDoubleArray.
     */
    public DataPipe getDataCaster(IDataInfo inputDataInfo) {
        if(inputDataInfo instanceof DataInfoDoubleArray) {
            return null;
        }
        return new CastToDoubleArray();
    }

    /**
     * Constructs the Data object used by this class.
     */
    private void setupData() {
        DataInfoDoubleArray independentInfo = new DataInfoDoubleArray(binnedDataInfo.getLabel(),binnedDataInfo.getDimension(),new int[]{histogram.getNBins()});
        data = new DataFunction(new int[]{histogram.getNBins()}, histogram.getHistogram());
        xDataSource = new DataSourceIndependentSimple(histogram.xValues(), independentInfo);
        dataInfo = new DataInfoFunction(binnedDataInfo.getLabel()+" Histogram", Null.DIMENSION, xDataSource);
        dataInfo.addTags(binnedDataInfo.getTags());
        dataInfo.addTag(getTag());
    }
    
    /**
     * @return the number of bins in the histogram.
     */
    public int getNBins() {
        return nBins;
    }

    /**
     * Sets the number of bins in each histogram. Calls setNBins method of
     * the histogram, which will discard data or modify itself
     * depending on how its own implementation.
     */
    public void setNBins(int nBins) {
        this.nBins = nBins;
        histogram.setNBins(nBins);
        setupData();
        if(dataSink != null) {
            dataSink.putDataInfo(dataInfo);
        }
    }

    public Histogram getHistogram() {
        return histogram;
    }

    /**
     * Returns the DataInfo for the output Data.
     */
    public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    private static final long serialVersionUID = 1L;
    protected Histogram histogram;
    protected DataSourceIndependentSimple xDataSource;
    private DataFunction data;
    private IDataInfo binnedDataInfo;
    protected int nData;
    private Histogram.Factory histogramFactory;
    private int nBins;
}
