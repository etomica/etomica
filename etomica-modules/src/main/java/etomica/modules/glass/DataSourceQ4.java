/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Time;

import java.util.Arrays;

public class DataSourceQ4 implements IDataSource, ConfigurationStorage.ConfigurationStorageListener, DataSourceIndependent {

    protected final ConfigurationStorage configStorage;
    protected DataDoubleArray tData;
    protected DataDoubleArray.DataInfoDoubleArray tDataInfo;
    protected DataFunction data;
    protected DataFunction.DataInfoFunction dataInfo;
    protected DataFunction chi4Data;
    protected DataFunction.DataInfoFunction chi4DataInfo;
    protected double[] Q4sum, Q4sum2;
    protected final DataTag tTag, tag, chi4Tag;

    protected long[] nSamples;
    protected int log2StepMin;
    protected double maxDr2 = 0.2 * 0.2;
    protected final Vector[] r;
    protected final int[] mobileAtoms;
    protected final Vector dr;

    public DataSourceQ4(ConfigurationStorage configStorage, int log2StepMin) {
        this.configStorage = configStorage;

        this.log2StepMin = log2StepMin;
        Q4sum2 = Q4sum = new double[0];
        nSamples = new long[0];
        tag = new DataTag();
        chi4Tag = new DataTag();
        tTag = new DataTag();
        int numAtoms = configStorage.getBox().getLeafList().size();
        Space space = configStorage.getBox().getSpace();
        r = space.makeVectorArray(numAtoms);
        dr = space.makeVector();
        mobileAtoms = new int[numAtoms];
        reallocate(0);
    }

    public double getMaxDr() {
        return Math.sqrt(maxDr2);
    }

    public void setMaxDr(double newMaxDr) {
        maxDr2 = newMaxDr * newMaxDr;
        zeroData();
    }

    public void zeroData() {
        reallocate(0);
    }

    public void reallocate(int n) {
        Q4sum = Arrays.copyOf(Q4sum, n);
        Q4sum2 = Arrays.copyOf(Q4sum2, n);
        nSamples = Arrays.copyOf(nSamples, n);
        data = new DataFunction(new int[]{n});
        chi4Data = new DataFunction(new int[]{n});
        tData = new DataDoubleArray(new int[]{n});
        tDataInfo = new DataDoubleArray.DataInfoDoubleArray("t", Time.DIMENSION, new int[]{n});
        tDataInfo.addTag(tTag);
        dataInfo = new DataFunction.DataInfoFunction("Q4", Null.DIMENSION, this);
        dataInfo.addTag(tag);
        chi4DataInfo = new DataFunction.DataInfoFunction("Q4^2", Null.DIMENSION, this);
        chi4DataInfo.addTag(chi4Tag);

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
    public void newConfigruation() {
        long step = configStorage.getSavedSteps()[0];
        Vector[] positions = configStorage.getSavedConfig(0);
        int numAtoms = positions.length;
        Boundary boundary = configStorage.getBox().getBoundary();
        for (int i = 0; i < configStorage.getLastConfigIndex(); i++) {
            int x = Math.max(i, log2StepMin);
            if (step % (1L << x) == 0) {
                if (i >= Q4sum.length) reallocate(i + 1); // reallocates if needed
                Vector[] iPositions = configStorage.getSavedConfig(i + 1);
                int nMobile = 0;
                int nReplaced = 0;
                for (int j = 0; j < numAtoms; j++) {
                    double r2 = positions[j].Mv1Squared(iPositions[j]);
                    if (r2 > maxDr2) {
                        mobileAtoms[nMobile] = j;
                        nMobile++;
                    }
                }
                // now look at all mobile to see if they've replaced
                for (int j = 0; j < nMobile; j++) {
                    int jj = mobileAtoms[j];
                    for (int k = 0; k < nMobile; k++) {
                        int kk = mobileAtoms[k];
                        dr.Ev1Mv2(positions[jj], iPositions[kk]);
                        boundary.nearestImage(dr);
                        double r2 = dr.squared();
                        if (r2 < maxDr2) {
                            nReplaced++;
                            // don't need to look for j replacing another
                            break;
                        }
                    }
                }
                double Q4 = ((numAtoms - nMobile) + nReplaced) / (double) numAtoms;
                Q4sum[i] += Q4;
                Q4sum2[i] += Q4 * Q4;
                nSamples[i]++;
            }
        }
    }


    @Override
    public IData getData() {
        if (configStorage.getLastConfigIndex() < 1) return data;
        double[] y = data.getData();
        for (int i = 0; i < Q4sum.length; i++) {
            y[i] = Q4sum[i] / (nSamples[i]);
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

    public MeterChi4 makeChi4Meter() {
        return new MeterChi4();
    }

    public class MeterChi4 implements IDataSource {

        @Override
        public IData getData() {
            if (configStorage.getLastConfigIndex() < 1) return chi4Data;
            double[] y = chi4Data.getData();
            for (int i = 0; i < Q4sum2.length; i++) {
                double avg = Q4sum[i] / nSamples[i];
                y[i] = Q4sum2[i] / nSamples[i] - avg * avg;
            }
            return chi4Data;
        }

        @Override
        public DataTag getTag() {
            return chi4Tag;
        }

        @Override
        public IDataInfo getDataInfo() {
            return chi4DataInfo;
        }
    }
}
