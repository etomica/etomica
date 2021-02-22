/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.data.*;
import etomica.data.history.HistoryScrolling;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.units.dimensions.Null;

import java.util.Arrays;

public class DataSourcePMSDHistory implements IDataSink, IDataSource, DataSourceIndependent, DataSourceMSD.MSDSink {

    protected DataDoubleArray xData;
    protected DataDoubleArray.DataInfoDoubleArray xDataInfo;
    protected DataFunction data;
    protected DataFunction.DataInfoFunction dataInfo;
    protected final DataTag xTag, tag;
    protected double[] lastSampleMSD, lastSampleX;
    protected double[] blockSumX;
    protected long[] nBlockSamplesX;
    protected long[] lastStepMSD, lastStepX;
    protected final TimeSource timeSource;
    protected long step0;
    protected boolean enabled;
    protected final int xlog2Interval;
    protected HistoryScrolling[] histories;
    protected int interval;

    public DataSourcePMSDHistory(TimeSource timeSource) {
        this(timeSource, 0);
    }

    public DataSourcePMSDHistory(TimeSource timeSource, int xlog2Interval) {
        this.timeSource = timeSource;
        this.xlog2Interval = xlog2Interval;
        lastSampleMSD = lastSampleX = new double[0];
        nBlockSamplesX = new long[60];
        lastStepX = new long[60];
        blockSumX = new double[60];
        histories = new HistoryScrolling[0];
        lastStepMSD = new long[0];
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
        }
        reallocate(0);
        setInterval(interval);
    }

    public void setInterval(int i) {
        interval = i;
        if (histories.length > interval && enabled) {
            double[] y = histories[i].getHistory();
            data = new DataFunction(new int[]{y.length}, y);
            double[] x = histories[i].getXValues();
            xData = new DataDoubleArray(new int[]{x.length}, x);
            xDataInfo = new DataDoubleArray.DataInfoDoubleArray("X", Null.DIMENSION, new int[]{x.length});
            dataInfo = new DataFunction.DataInfoFunction("MSD", Null.DIMENSION, this);
        } else {
            data = new DataFunction(new int[]{0});
            xData = new DataDoubleArray(new int[]{0});
            xDataInfo = new DataDoubleArray.DataInfoDoubleArray("X", Null.DIMENSION, new int[]{0});
            dataInfo = new DataFunction.DataInfoFunction("MSD", Null.DIMENSION, this);
        }
        xDataInfo.addTag(xTag);
        dataInfo.addTag(tag);
    }

    public void reallocate(int n) {
        if (n == lastSampleMSD.length && data != null) return;
        lastStepMSD = Arrays.copyOf(lastStepMSD, n);
        lastSampleX = Arrays.copyOf(lastSampleX, n);
        lastSampleMSD = Arrays.copyOf(lastSampleMSD, n);
        int nold = histories.length;
        histories = Arrays.copyOf(histories, n);
        for (int i = nold; i < n; i++) {
            histories[i] = new HistoryScrolling(500);
        }
        if (n == interval + 1) setInterval(interval);
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
        if (interval < histories.length) {
            histories[interval].getHistory();
            histories[interval].getXValues();
        }
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
                if (lastSampleX.length <= i) reallocate(i + 1);
                lastSampleX[i] = blockSumX[i] / nBlockSamplesX[i];
                lastStepX[i] = step;
                blockSumX[i] = 0;
                nBlockSamplesX[i] = 0;
                if (lastStepMSD[i] == step) {
                    processSample(i);
                }
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

    public void putMSD(int log2interval, long step, double msd) {
        if (log2interval < xlog2Interval) return;
        log2interval -= xlog2Interval;
        if (lastSampleMSD.length <= log2interval) reallocate(log2interval + 1);
        lastSampleMSD[log2interval] = msd;
        if (lastStepX[log2interval] == step) {
            processSample(log2interval);
        }
    }

    private void processSample(int log2interval) {
        histories[log2interval].addValue(lastSampleX[log2interval], lastSampleMSD[log2interval]);
    }
}
