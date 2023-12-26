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
 * Computes correlation of two input data streams, expected to each be 2 values.
 * Correlation is computed as though the two input data are vectors.
 *
 * When used to compute correlation between a dynamic (mobility) structure factor
 * and a static (density) structure factor from a DataSinkBlockAveragerSFac, the
 * effect is that the dynamic sfac is computed between slightly different configs than
 * those for the static sfac.  For an interval of 4, the static sfac will be averaged
 * over configs from steps 1, 2, 3 and 4.  The mobility will use the displacement from
 * steps 0 to 4.
 */
public class DataSourceCorrelation implements DataSourceIndependent {

    protected final ConfigurationStorage configStorage;
    protected DataDoubleArray tData;
    protected DataDoubleArray.DataInfoDoubleArray tDataInfo;
    protected final DataFunction[] data;
    protected DataFunction.DataInfoFunction[] dataInfo;
    protected double[][] xySum, x2Sum, y2Sum;
    protected final DataTag[] tag;
    protected final DataTag tTag;
    protected long[] nSamples;
    protected double[][][] lastSampleX, lastSampleY;
    protected long[] lastStepX, lastStepY;
    protected final int nData;

    public DataSourceCorrelation(ConfigurationStorage configStorage, int nData) {
        this.nData = nData;
        this.configStorage = configStorage;
        lastSampleX = new double[nData][0][0];
        lastSampleY = new double[nData][0][0];
        xySum = new double[nData][0];
        x2Sum = new double[nData][0];
        y2Sum = new double[nData][0];
        lastStepX = lastStepY = nSamples = new long[0];
        tag = new DataTag[nData];
        for (int i = 0; i < nData; i++) {
            tag[i] = new DataTag();
        }
        tTag = new DataTag();
        data = new DataFunction[nData];
        dataInfo = new DataFunction.DataInfoFunction[nData];
        reallocate(0);
    }

    public void reset() {
        reallocate(0);
    }

    public void reallocate(int n) {
        int nOld = nSamples.length;
        nSamples = Arrays.copyOf(nSamples, n);
        lastStepX = Arrays.copyOf(lastStepX, n);
        lastStepY = Arrays.copyOf(lastStepY, n);
        tData = new DataDoubleArray(new int[]{n});
        tDataInfo = new DataDoubleArray.DataInfoDoubleArray("t", Time.DIMENSION, new int[]{n});
        tDataInfo.addTag(tTag);
        for (int j = 0; j < nData; j++) {
            x2Sum[j] = Arrays.copyOf(x2Sum[j], n);
            y2Sum[j] = Arrays.copyOf(y2Sum[j], n);
            xySum[j] = Arrays.copyOf(xySum[j], n);
            lastSampleX[j] = Arrays.copyOf(lastSampleX[j], n);
            lastSampleY[j] = Arrays.copyOf(lastSampleY[j], n);
            for (int i = nOld; i < n; i++) {
                lastSampleX[j][i] = new double[2];
                lastSampleY[j][i] = new double[2];
            }
            data[j] = new DataFunction(new int[]{n});
            dataInfo[j] = new DataFunction.DataInfoFunction("X-Y cor", Null.DIMENSION, this);
            dataInfo[j].addTag(tag[j]);
        }
        setTimeData();
    }

    protected void setTimeData() {
        double[] t = tData.getData();
        if (t.length > 0) {
            double dt = configStorage.getSavedTimes()[0] - configStorage.getSavedTimes()[1];
            if (dt < 0) {
                Arrays.fill(t, Double.NaN);
            } else {
                for (int i = 0; i < t.length; i++) {
                    t[i] = dt * (1L << i);
                }
            }
        }
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

    private static double mydot(double[] x1, double[] x2) {
        return x1[0] * x2[0] + x1[1] * x2[1];
    }

    public void putData(int idx, int interval, double[][] xyData) {
        if (Double.isNaN(xyData[0][0])) return;
        long step = configStorage.getSavedSteps()[0];
        if (interval >= lastSampleX[0].length) reallocate(interval + 1);
        else if (Double.isNaN(tData.getValue(0))) setTimeData();
        for (int i = 0; i < nData; i++) {
            if (idx == 0) {
                lastSampleX[i][interval][0] = xyData[i][0];
                lastSampleX[i][interval][1] = xyData[i][1];
                if (i == 0) lastStepX[interval] = step;
            } else {
                lastSampleY[i][interval][0] = xyData[i][0];
                lastSampleY[i][interval][1] = xyData[i][1];
                if (i == 0) lastStepY[interval] = step;
            }
        }
        if (lastStepX[interval] == lastStepY[interval]) {
            for (int i = 0; i < nData; i++) {
                xySum[i][interval] += mydot(lastSampleX[i][interval], lastSampleY[i][interval]);
                x2Sum[i][interval] += mydot(lastSampleX[i][interval], lastSampleX[i][interval]);
                y2Sum[i][interval] += mydot(lastSampleY[i][interval], lastSampleY[i][interval]);
            }
            nSamples[interval]++;
            lastStepX[interval] = lastStepY[interval] = -1;
        }
    }

    public Receiver makeReceiver(int dataIdx) {
        return new Receiver(dataIdx);
    }

    public Meter makeMeter(int idx) {
        return new Meter(idx);
    }

    public class Receiver implements DataSinkBlockAveragerSFac.Sink {
        protected final int idx;

        public Receiver(int idx) {
            this.idx = idx;
        }

        public void putData(int interval, double[][] xyData) {
            DataSourceCorrelation.this.putData(idx, interval, xyData);
        }
    }

    public class Meter implements IDataSource {

        protected final int idx;

        public Meter(int idx) {
            this.idx = idx;
        }

        @Override
        public IData getData() {
            double[] y = data[idx].getData();
            for (int i = 0; i < y.length - 1; i++) {
                if (nSamples[i] < 2) {
                    y[i] = Double.NaN;
                    continue;
                }
                y[i] = xySum[idx][i] / Math.sqrt(x2Sum[idx][i] * y2Sum[idx][i]);
            }
            if (y.length > 0) y[y.length - 1] = Double.NaN;
            return data[idx];
        }

        @Override
        public DataTag getTag() {
            return tag[idx];
        }

        @Override
        public IDataInfo getDataInfo() {
            return dataInfo[idx];
        }
    }
}
