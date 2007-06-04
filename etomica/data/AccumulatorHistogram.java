package etomica.data;

import etomica.data.types.DataFunction;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.units.Null;
import etomica.util.Histogram;
import etomica.util.HistogramCollapsing;

/**
 * Accumulator that keeps histogram of data.
 * <p>
 * Input Data must implement DataArithmetic.
 */
public class AccumulatorHistogram extends DataAccumulator {

    /**
     * Creates instance using HistogramSimple factory and specifying histograms
     * having 100 bins.
     */
    public AccumulatorHistogram() {
        this(new HistogramCollapsing());
    }

    /**
     * Creates instance using given histogram factory with default nBins of 100.
     */
    public AccumulatorHistogram(Histogram histogram) {
        this(histogram, 100);
    }

    /**
     * Creates instance using the given histogram factory making histograms having
     * the given number of bins.
     */
    public AccumulatorHistogram(Histogram histogram, int nBins) {
        this.nBins = nBins;
        this.histogram = histogram;
    }

    /**
     * Adds each value in the given Data to its own histogram.
     */
    protected void addData(Data inputData) {
        histogram.addValue(inputData.getValue(0));
    }

    /**
     * Returns the set of histograms.
     */
    public Data getData() {
        // check to see if current data functions are the right length
        // histogram might change the number of bins on its own.
        // covertly call getHistogram() here.  We have to call it anyway
        // so that the array data gets updated.
        if (data.getData() != histogram.getHistogram() ||
                xDataSource.getIndependentData(0).getData() != histogram.xValues()) {
            setupData();
        }
        
        return data;
    }
    
    /**
     * Sets up data and histograms, discarding any previous results.
     */
    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        binnedDataInfo = inputDataInfo;
        if (inputDataInfo.getLength() != 1) {
            throw new IllegalArgumentException("AccumulatorHistogram can only handle single data");
        }
        setupData();
        return dataInfo;
    }
    
    /**
     * Returns null.  AccumulatorHistory can take an type of Data.
     */
    public DataPipe getDataCaster(IDataInfo inputDataInfo) {
        return null;
    }

    /**
     * Constructs the Data objects used by this class.
     */
    private void setupData() {
        DataInfoDoubleArray independentInfo = new DataInfoDoubleArray(binnedDataInfo.getLabel(),binnedDataInfo.getDimension(),new int[]{histogram.getNBins()});
        data = new DataFunction(new int[]{histogram.getNBins()}, histogram.getHistogram());
        xDataSource = new DataSourceIndependentSimple(histogram.xValues(), independentInfo);
        dataInfo = new DataInfoFunction(binnedDataInfo.getLabel()+" Histogram", Null.DIMENSION, xDataSource);
        dataInfo.addTags(binnedDataInfo.getTags());
        dataInfo.addTag(getTag());
        if (dataSink != null) {
            dataSink.putDataInfo(dataInfo);
        }
    }
    
    /**
     * @return the number of bins in each histogram.
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
        histogram.setNBins(nBins);
        setupData();
        if(dataSink != null) {
            dataSink.putDataInfo(dataInfo);
        }
    }

    /**
     * Zeros histograms, discarding any previous contributions. 
     */
    public void reset() {
        histogram.reset();
        setupData();
    }
    
    public Histogram getHistograms() {
        return histogram;
    }

    public void setHistogram(Histogram newHistogram) {
        histogram = newHistogram;
        reset();
    }

    /**
     * Returns the DataInfo for the output Data.
     */
    public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    private static final long serialVersionUID = 2L;
    protected Histogram histogram;
    protected DataSourceIndependentSimple xDataSource;
    private DataFunction data;
    private IDataInfo binnedDataInfo;
    private int nBins;
}
