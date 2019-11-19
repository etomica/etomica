/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.space.Vector;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Time;

import java.util.Arrays;

/**
 * Computes the excess kurtosis (alpha2) for the distribution of displacements
 */
public class DataSourceAlpha2 implements IDataSource, ConfigurationStorage.ConfigurationStorageListener, DataSourceIndependent {

    protected final ConfigurationStorage configStorage;
    protected DataDoubleArray tData;
    protected DataDoubleArray.DataInfoDoubleArray tDataInfo;
    protected DataFunction data;
    protected DataFunction.DataInfoFunction dataInfo;
    protected double[] msdSum, m4dSum;
    protected final DataTag tTag, tag;
    protected long[] nSamples;

    public DataSourceAlpha2(ConfigurationStorage configStorage) {
        this.configStorage = configStorage;
        msdSum = new double[0];
        m4dSum = new double[0];
        nSamples = new long[0];
        tag = new DataTag();
        tTag = new DataTag();
        reset();
    }

    public void reset() {
        int n = configStorage.getLastConfigIndex();
        if (n + 1 == msdSum.length && data != null) return;
        if (n < 1) n = 0;
        else n--;
        msdSum = Arrays.copyOf(msdSum, n);
        m4dSum = Arrays.copyOf(m4dSum, n);
        nSamples = Arrays.copyOf(nSamples, n);
        data = new DataFunction(new int[]{n});
        tData = new DataDoubleArray(new int[]{n});
        tDataInfo = new DataDoubleArray.DataInfoDoubleArray("t", Time.DIMENSION, new int[]{n});
        dataInfo = new DataFunction.DataInfoFunction("alpha", Null.DIMENSION, this);
        tDataInfo.addTag(tTag);
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

    @Override
    public IData getData() {
        if (configStorage.getLastConfigIndex() < 1) return data;
        double[] y = data.getData();
        int nAtoms = configStorage.getSavedConfig(0).length;
        // (3/5) for 3D; (1/2) for 2D
        double fac = configStorage.getBox().getSpace().D() == 2 ? 0.5 : 0.6;
        for (int i = 0; i < msdSum.length; i++) {
            y[i] = fac * m4dSum[i] / (msdSum[i] * msdSum[i]) * (nAtoms * nSamples[i]) - 1;
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
        reset(); // reallocates if needed
        long step = configStorage.getSavedSteps()[0];
        Vector[] positions = configStorage.getSavedConfig(0);
        for (int i = 1; i < msdSum.length; i++) {
            if (step % (1L << (i - 1)) == 0) {
                Vector[] iPositions = configStorage.getSavedConfig(i);
                for (int j = 0; j < positions.length; j++) {
                    double d2 = positions[j].Mv1Squared(iPositions[j]);
                    msdSum[i - 1] += d2;
                    m4dSum[i - 1] += d2 * d2;
                }
                nSamples[i - 1]++;
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
}
