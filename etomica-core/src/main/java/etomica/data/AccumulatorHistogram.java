/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.histogram.Histogram;
import etomica.data.histogram.HistogramCollapsing;
import etomica.data.histogram.HistogramNotSoSimple;
import etomica.data.histogram.HistogramReweightedData;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.units.dimensions.Null;

/**
 * Accumulator that keeps histogram of data.
 * <p>
 * Input Data must implement DataArithmetic.
 */
public class AccumulatorHistogram extends DataAccumulator {

    protected Histogram histogram;
    protected DataSourceIndependentSimple xDataSource;
    private DataFunction data;
    private IDataInfo binnedDataInfo;
    private int nBins;
    
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
    protected boolean addData(IData inputData) {
    	if (histogram instanceof HistogramNotSoSimple) {
    		((HistogramNotSoSimple)histogram).addValue(inputData.getValue(0), inputData.getValue(1));
    	}
    	else if (histogram instanceof HistogramReweightedData) {
            ((HistogramReweightedData)histogram).addValue(inputData.getValue(0), inputData.getValue(1));
    	}
    	else {
    	    histogram.addValue(inputData.getValue(0));
    	}
    	return true;
    }

    /**
     * Returns the set of histograms.
     */
    public IData getData() {
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
        if (inputDataInfo.getLength() != 1 && inputDataInfo.getLength() != 2) {
            throw new IllegalArgumentException("AccumulatorHistogram can only handle single data");
        }
        setupData();
        return dataInfo;
    }

    /**
     * Constructs the Data objects used by this class.
     */
    private void setupData() {
        DataInfoDoubleArray independentInfo = new DataInfoDoubleArray(binnedDataInfo.getLabel(),binnedDataInfo.getDimension(),new int[]{histogram.getNBins()});
        data = new DataFunction(new int[]{histogram.getNBins()}, histogram.getHistogram());
        if (xDataSource != null) {
            xDataSource.update(histogram.xValues(), independentInfo);
        }
        else {
            xDataSource = new DataSourceIndependentSimple(histogram.xValues(), independentInfo);
        }
        dataInfo = new DataInfoFunction(binnedDataInfo.getLabel()+" Histogram", Null.DIMENSION, xDataSource);
        dataInfo.addTags(binnedDataInfo.getTags());
        dataInfo.addTag(tag);
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
}
