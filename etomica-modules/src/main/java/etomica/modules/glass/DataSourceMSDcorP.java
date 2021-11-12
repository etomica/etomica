/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.integrator.IntegratorMD;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Time;

import java.util.Arrays;

public class DataSourceMSDcorP implements IDataSink, IDataSource, DataSourceIndependent, DataSourceMSD.MSDSink {

    protected DataDoubleArray tData;
    protected DataDoubleArray.DataInfoDoubleArray tDataInfo;
    protected DataFunction data, covData;
    protected DataFunction.DataInfoFunction dataInfo, covDataInfo;
    protected double[] msdSum, xSum, xmsdSum, msd2Sum, x2Sum;
    protected final DataTag tTag, tag, covTag;
    protected long[] nSamples;
    protected double[] lastSampleMSD, lastSampleX;
    protected double[] blockSumX;
    protected long[] nBlockSamplesX;
    protected long[] lastStepMSD, lastStepX;
    protected final IntegratorMD integrator;
    protected long step0;
    protected boolean enabled;
    protected final int xlog2Interval;

    public DataSourceMSDcorP(IntegratorMD integrator) {
        this(integrator, 0);
    }

    public DataSourceMSDcorP(IntegratorMD integrator, int xlog2Interval) {
        this.integrator = integrator;
        this.xlog2Interval = xlog2Interval;
        lastSampleMSD = lastSampleX = xSum = xmsdSum = msdSum = msd2Sum = x2Sum = new double[0];
        nBlockSamplesX = new long[60];
        lastStepX = new long[60];
        blockSumX = new double[60];
        lastStepMSD = nSamples = new long[0];
        tag = new DataTag();
        tTag = new DataTag();
        covTag = new DataTag();
        resetStep0();
    }

    public void resetStep0() {
        step0 = -1;
        for (int i = 0; i < nBlockSamplesX.length; i++) {
            lastStepX[i] = nBlockSamplesX[i] = 0;
            blockSumX[i] = 0;
        }
        reallocate(0);
    }

    public void reallocate(int n) {
        if (n == msdSum.length && data != null) return;
        msdSum = Arrays.copyOf(msdSum, n);
        msd2Sum = Arrays.copyOf(msd2Sum, n);
        xSum = Arrays.copyOf(xSum, n);
        x2Sum = Arrays.copyOf(x2Sum, n);
        xmsdSum = Arrays.copyOf(xmsdSum, n);
        nSamples = Arrays.copyOf(nSamples, n);
        lastStepMSD = Arrays.copyOf(lastStepMSD, n);
        lastSampleX = Arrays.copyOf(lastSampleX, n);
        lastSampleMSD = Arrays.copyOf(lastSampleMSD, n);
        data = new DataFunction(new int[]{n});
        covData = new DataFunction(new int[]{n});
        tData = new DataDoubleArray(new int[]{n});
        tDataInfo = new DataDoubleArray.DataInfoDoubleArray("t", Time.DIMENSION, new int[]{n});
        tDataInfo.addTag(tTag);
        dataInfo = new DataFunction.DataInfoFunction("X-MSD cor", Null.DIMENSION, this);
        dataInfo.addTag(tag);
        covDataInfo = new DataFunction.DataInfoFunction("X-MSD cov", Null.DIMENSION, this);
        covDataInfo.addTag(covTag);
        double[] t = tData.getData();
        if (t.length > 0) {
            double dt = integrator.getTimeStep() * (1L << xlog2Interval);
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
        double[] y = data.getData();
        for (int i = 0; i < y.length - 1; i++) {
            double pAvg = xSum[i] / nSamples[i];
            double msdAvg = msdSum[i] / nSamples[i];
            double pmsdAvg = xmsdSum[i] / nSamples[i];
            double p2Avg = x2Sum[i] / nSamples[i];
            double pVar = p2Avg - pAvg * pAvg;
            double msd2Avg = msd2Sum[i] / nSamples[i];
            double msdVar = msd2Avg - msdAvg * msdAvg;
            y[i] = (pmsdAvg - pAvg * msdAvg) / Math.sqrt(msdVar * pVar);
        }
        if (y.length > 0) y[y.length - 1] = Double.NaN;
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
        if (step0 < 0) step0 = integrator.getStepCount() - (1L << xlog2Interval);
        long step = integrator.getStepCount() - step0;
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
        xSum[log2interval] += lastSampleX[log2interval];
        msdSum[log2interval] += lastSampleMSD[log2interval];
        xmsdSum[log2interval] += lastSampleX[log2interval] * lastSampleMSD[log2interval];
        x2Sum[log2interval] += lastSampleX[log2interval] * lastSampleX[log2interval];
        msd2Sum[log2interval] += lastSampleMSD[log2interval] * lastSampleMSD[log2interval];
        lastStepMSD[log2interval] = lastStepX[log2interval] = -1;
        nSamples[log2interval]++;
    }

    public DataSourceMSDcovP makeCov() {
        return new DataSourceMSDcovP();
    }

    /**
     * Returns the covariance of the MSD with P, divided by t
     */
    public class DataSourceMSDcovP implements IDataSource {

        @Override
        public IDataInfo getDataInfo() {
            return covDataInfo;
        }

        @Override
        public DataTag getTag() {
            return covTag;
        }

        @Override
        public IData getData() {
            double[] y = data.getData();
            for (int i = 0; i < y.length - 1; i++) {
                double pAvg = xSum[i] / nSamples[i];
                double msdAvg = msdSum[i] / nSamples[i];
                double pmsdAvg = xmsdSum[i] / nSamples[i];
                y[i] = pmsdAvg - pAvg * msdAvg / tData.getValue(i);
            }
            if (y.length > 0) y[y.length - 1] = Double.NaN;
            return data;
        }

    }
}
