/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.types.CastToDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.units.Null;
import etomica.math.DoubleRange;
import etomica.data.histogram.Histogram;
import etomica.data.histogram.HistogramSimple;

/**
 * DataProcessor that creates a histogram from each piece of data that comes in.
 * A single Histogram object is kept, but is reset for each new data.
 * <p>
 * Input Data must implement DataArithmetic.
 */
public class DataHistogram extends DataProcessor {

    /**
     * Creates instance using HistogramSimple and specifying histograms
     * having the given number of bins and the given range.
     */
    public DataHistogram(int nBins, DoubleRange range) {
        this(new HistogramSimple(nBins, range));
    }

    /**
     * Creates instance using the given histogram.
     */
    public DataHistogram(Histogram histogram) {
        this.histogram = histogram;
    }

    /**
     * Adds each value in the given Data to a single histogram.
     */
    protected IData processData(IData inputData) {
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
    protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
        binnedDataInfo = inputDataInfo;
        nData = inputDataInfo.getLength();
        setupData();
        return dataInfo;
    }
    
    /**
     * Returns caster that ensures accumulator will receive a DataDoubleArray.
     */
    public DataPipe getDataCaster(IEtomicaDataInfo inputDataInfo) {
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
        dataInfo.addTag(tag);
    }

    public Histogram getHistogram() {
        return histogram;
    }

    /**
     * Returns the DataInfo for the output Data.
     */
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    
    private static final long serialVersionUID = 1L;
    protected Histogram histogram;
    protected DataSourceIndependentSimple xDataSource;
    private DataFunction data;
    private IEtomicaDataInfo binnedDataInfo;
    protected int nData;
}
