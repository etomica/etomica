/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Time;

import java.util.Arrays;

/**
 * Computes the excess kurtosis (alpha2) for the distribution of displacements
 */
public class DataSourceF implements IDataSource, ConfigurationStorage.ConfigurationStorageListener, DataSourceIndependent {

    protected final ConfigurationStorage configStorage;
    protected DataDoubleArray tData;
    protected DataDoubleArray.DataInfoDoubleArray tDataInfo;
    protected DataFunction data;
    protected DataFunction.DataInfoFunction dataInfo;
    protected double[] fSum;
    protected final DataTag tTag, tag;
    protected long[] nSamples;
    protected final Vector dr, q;
    protected AtomType type;
    protected double strucFac =0;
    protected double[] cSum;
    protected double[] sSum;


    public DataSourceF(ConfigurationStorage configStorage) {
        this.configStorage = configStorage;
        Space space = configStorage.getBox().getSpace();
        fSum = new double[0];
        nSamples = new long[0];
        tag = new DataTag();
        tTag = new DataTag();
        dr = space.makeVector();
        q = space.makeVector();
        q.setX(0,7.0);
        cSum = new double[0];
        sSum = new double[0];
        reset();
    }

    public void reset() {
        int n = configStorage.getLastConfigIndex();
        if (n == fSum.length && data != null) return;
        if (n < 1) n = 0;
        if(n==0){
            strucFac = 0;
        }
        fSum = Arrays.copyOf(fSum, n);
        nSamples = Arrays.copyOf(nSamples, n);
        data = new DataFunction(new int[]{n});
        tData = new DataDoubleArray(new int[]{n});
        tDataInfo = new DataDoubleArray.DataInfoDoubleArray("t", Time.DIMENSION, new int[]{n});
        dataInfo = new DataFunction.DataInfoFunction("F(t)", Null.DIMENSION, this);
        dataInfo.addTag(tag);
        cSum = Arrays.copyOf(cSum, n+1);
        sSum = Arrays.copyOf(sSum, n+1);

        if(n != 0){
            cSum[n] = cSum[n-1];
            sSum[n] = sSum[n-1];
        }

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

        for (int i = 0; i < fSum.length; i++) {
            y[i] = fSum[i] / (strucFac/nSamples[0])/nSamples[i];
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
        double c0Sum =0 , s0Sum = 0;
        for (int j = 0; j < positions.length; j++) {
            double qdotr0 = q.dot(positions[j]);
            c0Sum += Math.cos(qdotr0);
            s0Sum += Math.sin(qdotr0);
        }
        strucFac += c0Sum*c0Sum + s0Sum*s0Sum;

        for (int i = 1, d = 1; d < step && (step-1) % d == 0; i++, d *= 2) {
            cSum[i] = cSum[0];
            sSum[i] = sSum[0];
        }

        cSum[0] = c0Sum;
        sSum[0] = s0Sum;

        for (int i = 1, d = 1; d < step + 1 && step % d == 0; i++, d *= 2) {
            fSum[i-1] += c0Sum*cSum[i-1] + s0Sum*sSum[i-1];
            nSamples[i-1]++;
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