/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.data.*;
import etomica.data.histogram.Histogram;
import etomica.data.histogram.HistogramCollapsing;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.units.dimensions.Null;

public class DataSourceHisogram implements IDataSink, IDataSource, DataSourceIndependent {

    protected DataDoubleArray xData;
    protected DataDoubleArray.DataInfoDoubleArray xDataInfo;
    protected DataFunction data;
    protected DataFunction.DataInfoFunction dataInfo;
    protected final DataTag xTag, tag;
    protected double[] blockSumX;
    protected long[] nBlockSamplesX;
    protected long[] lastStepX;
    protected final TimeSource timeSource;
    protected long step0;
    protected boolean enabled;
    protected final int xlog2Interval;
    protected int interval;
    protected Histogram[] histograms;

    public DataSourceHisogram(TimeSource timeSource) {
        this(timeSource, 0);
    }

    public DataSourceHisogram(TimeSource timeSource, int xlog2Interval) {
        this.timeSource = timeSource;
        this.xlog2Interval = xlog2Interval;
        histograms = new Histogram[60];
        for (int i = 0; i < histograms.length; i++) {
            histograms[i] = new HistogramCollapsing(100);
        }
        nBlockSamplesX = new long[60];
        lastStepX = new long[60];
        blockSumX = new double[60];
        tag = new DataTag();
        xTag = new DataTag();
        resetStep0();
    }

    public void resetStep0() {
        step0 = timeSource.getStepCount();
        step0 -= step0 % (1L << xlog2Interval);
        for (int i = 0; i < nBlockSamplesX.length; i++) {
            lastStepX[i] = nBlockSamplesX[i] = 0;
            blockSumX[i] = 0;
            histograms[i].reset();
        }
        setInterval(interval);
    }

    public void setInterval(int i) {
        interval = i;
        double[] y = histograms[i].getHistogram();
        data = new DataFunction(new int[]{y.length}, y);
        double[] x = histograms[i].xValues();
        xData = new DataDoubleArray(new int[]{x.length}, x);
        xDataInfo = new DataDoubleArray.DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{x.length});
        xDataInfo.addTag(xTag);
        dataInfo = new DataFunction.DataInfoFunction("stuff histogram", Null.DIMENSION, this);
        dataInfo.addTag(tag);
    }

    public int getInterval() {
        return interval;
    }

    public void setEnabled(boolean isEnabled) {
        enabled = isEnabled;
        resetStep0();
    }

    public boolean getEnabled() {
        return enabled;
    }

    @Override
    public IData getData() {
        histograms[interval].getHistogram();
        return data;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    @Override
    public void putData(IData inputData) {
        if (!enabled) return;
        double x = inputData.getValue(0);
        long step = timeSource.getStepCount() - step0;
        for (int i = 0; i < blockSumX.length; i++) {
            blockSumX[i] += x;
            nBlockSamplesX[i]++;
            if (step % (1L << (xlog2Interval + i)) == 0) {
                double sample = blockSumX[i] / nBlockSamplesX[i];
                histograms[i].addValue(sample);
                lastStepX[i] = step;
                blockSumX[i] = 0;
                nBlockSamplesX[i] = 0;
            }
        }
    }

    @Override
    public void putDataInfo(IDataInfo inputDataInfo) {
    }

    @Override
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    @Override
    public DataDoubleArray getIndependentData(int i) {
        return xData;
    }

    @Override
    public DataDoubleArray.DataInfoDoubleArray getIndependentDataInfo(int i) {
        return xDataInfo;
    }

    @Override
    public int getIndependentArrayDimension() {
        return 1;
    }

    @Override
    public DataTag getIndependentTag() {
        return xTag;
    }
}
