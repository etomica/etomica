/*
 * History
 * Created on Aug 4, 2004 by kofke
 */
package etomica.data;

import etomica.data.types.CastToDoubleArray;
import etomica.data.types.DataArithmetic;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataGroup;
import etomica.units.Dimension;
import etomica.util.Histogram;
import etomica.util.HistogramSimple;

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
        this(HistogramSimple.FACTORY);
    }

    /**
     * Creates instance using given histogram factory with default nBins of 100.
     */
    public AccumulatorHistogram(Histogram.Factory factory) {
        this(factory, 100);
    }

    /**
     * Creates instance using the given histogram factory making histograms having
     * the given number of bins.
     */
    public AccumulatorHistogram(Histogram.Factory factory, int nBins) {
        this.nBins = nBins;
        histogramFactory = factory;
    }

    /**
     * Adds each value in the given Data to its own histogram.
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
        // check to see if current data functions are the right length
        // histogram might change the number of bins on its own.
        boolean success = true;
        for (int i=0; i<nData; i++) {
            DataFunction dataFunction = (DataFunction)data.getData(i);
            if (dataFunction.getLength() != histogram[i].getNBins()) {
                success = false;
            }
        }
        if (!success) {
            DataFunction[] dataFunctions = new DataFunction[nData];
            // attempt to re-use old DataFunctions
            for (int i=0; i<nData; i++) {
                DataFunction dataFunction = (DataFunction)data.getData(i);
                if (dataFunction.getLength() == histogram[i].getNBins()) {
                    dataFunctions[i] = dataFunction;
                }
                else {
                    dataFunctions[i] = new DataFunction(binnedDataInfo.getLabel()+" Histogram",Dimension.NULL,
                            new DataDoubleArray[]{new DataDoubleArray(binnedDataInfo.getLabel(), 
                                    binnedDataInfo.getDimension(), histogram[i].getNBins())});
                }
            }
            // creating a new data instance might confuse downstream data sinks, but
            // we have little choice and they should deal.
            data = new DataGroup("Histogram",dataFunctions);
        }   
            
        for (int i = 0; i < nData; i++) {
            DataFunction dataFunction = (DataFunction)data.getData(i);
            dataFunction.getYData().E(histogram[i].getHistogram());
            dataFunction.getXData(0).E(histogram[i].xValues());
        }
        return data;
    }
    
    /**
     * Sets up data and histograms, discarding any previous results.
     */
    protected DataInfo processDataInfo(DataInfo inputDataInfo) {
        binnedDataInfo = inputDataInfo;
        nData = ((DataDoubleArray.Factory)inputDataInfo.getDataFactory()).getArrayLength();
        histogram = new Histogram[nData];
        for (int i = 0; i < nData; i++) {
            histogram[i] = histogramFactory.makeHistogram(nBins);
        }
        setupData();
        return data.getDataInfo();
    }
    
    /**
     * Returns caster that ensures accumulator will receive a DataDoubleArray.
     */
    public DataProcessor getDataCaster(DataInfo inputDataInfo) {
        if(inputDataInfo.getDataClass() == DataDoubleArray.class) {
            return null;
        }
        return new CastToDoubleArray();
    }

    /**
     * Constructs the Data objects used by this class.
     */
    private void setupData() {
        DataFunction[] dataFunctions = new DataFunction[nData];
        for (int i=0; i<nData; i++) {
            dataFunctions[i] = new DataFunction(binnedDataInfo.getLabel()+" Histogram",Dimension.NULL,new DataDoubleArray[]{new DataDoubleArray(binnedDataInfo.getLabel(), binnedDataInfo.getDimension(), histogram[i].getNBins())});
        }
        data = new DataGroup("Histogram",dataFunctions);
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
        for (int i = 0; i < nData; i++) {
            histogram[i].setNBins(nBins);
        }
        setupData();
        if(dataSink != null) {
            dataSink.putDataInfo(data.getDataInfo());
        }
    }

    /**
     * Zeros histograms, discarding any previous contributions. 
     */
    public void reset() {
        for (int i = 0; i < nData; i++)
            histogram[i].reset();
    }

    /**
     * Returns the DataInfo for the output Data.
     */
    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }
    
    Histogram[] histogram = new Histogram[0];
    private DataGroup data;
    private DataInfo binnedDataInfo;
    int nData;
    private Histogram.Factory histogramFactory;
    private int nBins;

}
