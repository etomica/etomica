/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Time;

import java.util.Arrays;

/**
 * Measures correlation of blocks of increasing size.
 */
public class DataSourceBlockAvgCor implements IDataSink, IDataSource, DataSourceIndependent {

    protected DataDoubleArray tData;
    protected DataDoubleArray.DataInfoDoubleArray tDataInfo;
    protected DataFunction corPData;
    protected DataFunction.DataInfoFunction corPDataInfo;
    protected double[] xSum, x2Sum;
    protected final DataTag tTag, corPTag;
    protected long[] nSamplesCorX;
    protected double[] lastSampleX, firstSampleX;
    protected double[] blockSumX;
    protected double[] corSumX;
    protected final TimeSource timeSource;
    protected long step0;
    protected boolean enabled;
    protected int minInterval;

    public DataSourceBlockAvgCor(TimeSource timeSource) {
        this.timeSource = timeSource;
        blockSumX = x2Sum = xSum = firstSampleX = corSumX = lastSampleX = new double[0];
        nSamplesCorX = new long[0];
        tTag = new DataTag();
        corPTag = new DataTag();
        resetStep0();
    }

    public void setMinInterval(int minInterval) {
        this.minInterval = minInterval;
    }

    public void resetStep0() {
        step0 = timeSource.getStepCount();
        reallocate(0);
    }

    protected void reallocate(int n) {
        if (n == xSum.length && corPData != null) return;
        blockSumX = Arrays.copyOf(blockSumX, n);
        xSum = Arrays.copyOf(xSum, n);
        x2Sum = Arrays.copyOf(x2Sum, n);
        nSamplesCorX = Arrays.copyOf(nSamplesCorX, n);
        lastSampleX = Arrays.copyOf(lastSampleX, n);
        firstSampleX = Arrays.copyOf(firstSampleX, n);
        corSumX = Arrays.copyOf(corSumX, n);
        corPData = new DataFunction(new int[]{n});
        tData = new DataDoubleArray(new int[]{n});
        tDataInfo = new DataDoubleArray.DataInfoDoubleArray("t", Time.DIMENSION, new int[]{n});
        tDataInfo.addTag(tTag);
        corPDataInfo = new DataFunction.DataInfoFunction("P cor", Null.DIMENSION, this);
        corPDataInfo.addTag(corPTag);
        double[] t = tData.getData();
        if (t.length > 0) {
            double dt = timeSource.getTimeStep();
            for (int i = 0; i < t.length; i++) {
                t[i] = dt * (1L << i);
            }
        }
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
        double[] y = corPData.getData();
        for (int i = 0; i < y.length - 1; i++) {
            if (nSamplesCorX[i] < 4) {
                y[i] = Double.NaN;
                continue;
            }
            double avg1 = (xSum[i] - lastSampleX[i]) / nSamplesCorX[i];
            double avg2 = (xSum[i] - firstSampleX[i]) / nSamplesCorX[i];
            double avg21 = (x2Sum[i] - lastSampleX[i] * lastSampleX[i]) / nSamplesCorX[i];
            double avg22 = (x2Sum[i] - firstSampleX[i] * firstSampleX[i]) / nSamplesCorX[i];
            double avgCor = corSumX[i] / nSamplesCorX[i];

            // <(x2 - xa2)(x1 - xa1)> / < (x - xa)^2 >
            // (<x2 x1> - <x2 xa1> - <x1 xa2> + xa1 xa2) / ( <x^2> - <x xa> - <xa x> + <xa>^2)
            // (<x2 x1> - xa1 xa2) / (stdev1 stdev2)
            y[i] = (avgCor - avg1 * avg2) / Math.sqrt((avg21 - avg1 * avg1) * (avg22 - avg2 * avg2));
        }
        if (y.length > 0) y[y.length - 1] = Double.NaN;
        return corPData;
    }

    @Override
    public DataTag getTag() {
        return corPTag;
    }

    @Override
    public void putData(IData inputData) {
        if (!enabled) return;
        double x = inputData.getValue(0);
        if (blockSumX.length == 0) step0 = timeSource.getStepCount() - (1L << minInterval);
        long step = timeSource.getStepCount() - step0;
        if (blockSumX.length == 0) reallocate(1);
        for (int i = minInterval; step >= 1L << i; i++) {
            x = processData(i, step, x);
            if (Double.isNaN(x)) break;
        }
    }

    protected double processData(int interval, long step, double blockAvg) {
        if (blockSumX.length <= interval) reallocate(interval + 1);
        if (step == (1L << interval)) {
            firstSampleX[interval] = blockAvg;
        } else {
            corSumX[interval] += blockAvg * lastSampleX[interval];
            nSamplesCorX[interval]++;
        }
        xSum[interval] += blockAvg;
        x2Sum[interval] += blockAvg * blockAvg;
        double rv = Double.NaN;
        if (step % (1L << (interval + 1)) == 0) rv = (lastSampleX[interval] + blockAvg) / 2.0;
        lastSampleX[interval] = blockAvg;
        return rv;
    }

    @Override
    public void putDataInfo(IDataInfo inputDataInfo) {
    }

    @Override
    public IDataInfo getDataInfo() {
        return corPDataInfo;
    }

    @Override
    public DataDoubleArray getIndependentData(int i) {
        return tData;
    }

    @Override
    public DataDoubleArray.DataInfoDoubleArray getIndependentDataInfo(int i) {
        return tDataInfo;
    }

    @Override
    public int getIndependentArrayDimension() {
        return 1;
    }

    @Override
    public DataTag getIndependentTag() {
        return tTag;
    }
}
