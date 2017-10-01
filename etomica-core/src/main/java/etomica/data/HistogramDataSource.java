/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.histogram.Histogram;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.units.dimensions.Null;

/**
 * DataSource that simply exposes a Histogram as a DataFunction.
 */
public class HistogramDataSource implements IDataSource {

    protected final Histogram histogram;
    protected DataSourceIndependentSimple xDataSource;
    private DataFunction data;
    private IEtomicaDataInfo dataInfo;
    protected final DataTag tag;

    public HistogramDataSource(Histogram histogram) {
        this.histogram = histogram;
        tag = new DataTag();
        setupData();
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
     * Constructs the Data objects used by this class.
     */
    private void setupData() {
        DataInfoDoubleArray independentInfo = new DataInfoDoubleArray("x",Null.DIMENSION,new int[]{histogram.getNBins()});
        data = new DataFunction(new int[]{histogram.getNBins()}, histogram.getHistogram());
        if (xDataSource != null) {
            xDataSource.update(histogram.xValues(), independentInfo);
        }
        else {
            xDataSource = new DataSourceIndependentSimple(histogram.xValues(), independentInfo);
        }
        dataInfo = new DataInfoFunction("histogram", Null.DIMENSION, xDataSource);
        dataInfo.addTag(tag);
    }

    public DataTag getTag() {
        return tag;
    }

    /**
     * Returns the DataInfo for the output Data.
     */
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
}
