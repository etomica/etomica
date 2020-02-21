/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.atom.AtomType;
import etomica.atom.IAtomList;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.space.Vector;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Time;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class DataSourceMSD implements IDataSource, ConfigurationStorage.ConfigurationStorageListener, DataSourceIndependent {

    protected final ConfigurationStorage configStorage;
    protected DataDoubleArray tData;
    protected DataDoubleArray.DataInfoDoubleArray tDataInfo;
    protected DataFunction data, errData;
    protected DataFunction.DataInfoFunction dataInfo;
    protected double[] msdSum, msd2Sum, msdSumBlock;
    protected final DataTag tTag, tag;
    protected long[] nSamples;
    protected final AtomType type;
    protected List<MSDSink> msdSinks;
    protected int minInterval = 3;

    public DataSourceMSD(ConfigurationStorage configStorage) {
        this(configStorage, null);
    }

    public DataSourceMSD(ConfigurationStorage configStorage, AtomType type) {
        this.configStorage = configStorage;
        this.type = type;
        msdSum = new double[0];
        msd2Sum = new double[0];
        msdSumBlock = new double[0];
        nSamples = new long[0];
        tag = new DataTag();
        tTag = new DataTag();
        msdSinks = new ArrayList<>();
        reallocate(0);
    }

    public void reset() {
        reallocate(0);
    }

    protected void reallocate(int n) {
        msdSum = Arrays.copyOf(msdSum, n);
        msd2Sum = Arrays.copyOf(msd2Sum, n);
        msdSumBlock = Arrays.copyOf(msdSumBlock, n);
        nSamples = Arrays.copyOf(nSamples, n);
        data = new DataFunction(new int[]{n});
        errData = new DataFunction(new int[]{n});
        tData = new DataDoubleArray(new int[]{n});
        tDataInfo = new DataDoubleArray.DataInfoDoubleArray("t", Time.DIMENSION, new int[]{n});
        tDataInfo.addTag(tTag);
        dataInfo = new DataFunction.DataInfoFunction("MSD", new CompoundDimension(new Dimension[]{Length.DIMENSION}, new double[]{2}), this);
        dataInfo.addTag(tag);
        double[] t = tData.getData();
        if (t.length > 0) {
            double[] savedTimes = configStorage.getSavedTimes();
            double dt = savedTimes[0] - savedTimes[1];
            for (int i = 0; i < t.length; i++) {
                t[i] = dt * (1L << i);
            }
        }
    }

    public void addMSDSink(MSDSink sink) {
        msdSinks.add(sink);
    }

    @Override
    public IData getData() {
        if (configStorage.getLastConfigIndex() < 1) return data;
        double[] y = data.getData();
        double[] yErr = errData.getData();
        for (int i = 0; i < msdSum.length; i++) {
            long M = nSamples[i];
            y[i] = msdSum[i] / M;
            yErr[i] = Math.sqrt((msd2Sum[i]/M - y[i]*y[i]) / (M - 1));
        }
        return data;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    @Override
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    @Override
    public void newConfigruation() {
        int blockSize = 1;
        long step = configStorage.getSavedSteps()[0];
        Vector[] positions = configStorage.getSavedConfig(0);
        IAtomList atoms = configStorage.getBox().getLeafList();
        for (int i = 0; i < configStorage.getLastConfigIndex(); i++) {
            int x = Math.max(i, minInterval);
            if (step % (1L << x) == 0) {
                if (i >= msdSum.length) reallocate(i + 1);
                Vector[] iPositions = configStorage.getSavedConfig(i + 1);
                double iSum = 0;
                int iSamples = 0;
                for (int j = 0; j < positions.length; j++) {
                    if (type != null && atoms.get(j).getType() != type) continue;
                    iSum += positions[j].Mv1Squared(iPositions[j]);
                    iSamples++;
                }
                double iAvg = iSum / iSamples;
                msdSumBlock[i] += iAvg;
                if (step % (blockSize * (1L << i)) == 0) {
                    double xb = msdSumBlock[i] / blockSize;
                    msdSum[i] += xb;
                    msd2Sum[i] += xb * xb;
                    nSamples[i]++;
                    msdSumBlock[i] = 0;
                }
                for (MSDSink s : msdSinks) {
                    s.putMSD(i, step, iAvg);
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

    public interface MSDSink {
        void putMSD(int log2interval, long step, double msd);
    }
}
