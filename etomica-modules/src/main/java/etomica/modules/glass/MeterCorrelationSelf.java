/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.glass;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Time;

import java.util.Arrays;

public class MeterCorrelationSelf implements ConfigurationStorage.ConfigurationStorageListener, IDataSource, DataSourceIndependent {

    public enum CorrelationType {TOTAL, MAGNITUDE, MAG_DOT}

    protected CorrelationType correlationType;
    protected final ConfigurationStorage configStorage;
    protected final DataTag tag, tTag;
    protected final Vector dr01, dr12;
    protected double[] corSum;
    protected long[] nSamples;
    protected DataFunction data;
    protected DataDoubleArray tData;
    protected DataDoubleArray.DataInfoDoubleArray tDataInfo;
    protected AtomType type;
    private IDataInfo dataInfo;
    protected double[] dr1Sum, dr2Sum, dr3Sum;
    protected int minInterval = 3;

    public MeterCorrelationSelf(ConfigurationStorage configStorage) {
        this(configStorage, CorrelationType.TOTAL);
    }

    public MeterCorrelationSelf(ConfigurationStorage configStorage, CorrelationType cType) {

        this.configStorage = configStorage;
        this.correlationType = cType;
        Space space = configStorage.getBox().getSpace();
        dr1Sum = dr2Sum = dr3Sum = corSum = new double[0];
        nSamples = new long[0];

        dr01 = space.makeVector();
        dr12 = space.makeVector();
        tag = new DataTag();
        tTag = new DataTag();
        reset();
    }

    public void reset() {
        int n = configStorage.getLastConfigIndex();
        if (n + 1 == corSum.length && data != null) return;
        if (n < 1) n = 0;
        else n--;
        corSum = Arrays.copyOf(corSum, n);
        dr1Sum = Arrays.copyOf(dr1Sum, n);
        dr2Sum = Arrays.copyOf(dr2Sum, n);
        dr3Sum = Arrays.copyOf(dr3Sum, n);
        nSamples = Arrays.copyOf(nSamples, n);
        data = new DataFunction(new int[]{n});
        tData = new DataDoubleArray(new int[]{n});
        tDataInfo = new DataDoubleArray.DataInfoDoubleArray("t", Time.DIMENSION, new int[]{n});
        dataInfo = new DataFunction.DataInfoFunction("cor", Null.DIMENSION, this);
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

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public void setAtomType(AtomType type) {
        this.type = type;
    }

    /**
     * Returns the RDF, averaged over the calls to actionPerformed since the
     * meter was reset or had some parameter changed (xMax or # of bins).
     */
    public IData getData() {
        if (configStorage.getLastConfigIndex() < 2) {
            data.E(Double.NaN);
            return data;
        }

        final double[] y = data.getData();
        for (int i = 0; i < y.length; i++) {
            double avg1 = dr1Sum[i] / (nSamples[i] * 2);
            double avg2 = dr2Sum[i] / (nSamples[i] * 2);
            double avg3 = dr3Sum[i] / (nSamples[i] * 2);
            if (correlationType == CorrelationType.TOTAL) {
                y[i] = (corSum[i] / nSamples[i]) / avg2;
            } else if (correlationType == CorrelationType.MAGNITUDE) {
                y[i] = (corSum[i] / nSamples[i] - avg1 * avg1) /
                        (avg2 - avg1 * avg1);

            } else if (correlationType == CorrelationType.MAG_DOT) {
                y[i] = (corSum[i] / (2 * nSamples[i])) / avg3;
            }
        }
        return data;
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

    @Override
    public void newConfigruation() {

        reset(); // reallocates if needed
        long step2 = configStorage.getSavedSteps()[0];
        Vector[] config2 = configStorage.getSavedConfig(0);
        IAtomList atoms = configStorage.getBox().getLeafList();

        // j corresponds to the 0 -> 2 interval, so we need to start
        // with j=1 (step=2 and we look at steps 0,1,2) and then store in j-1
        for (int j = 1; j < configStorage.getLastConfigIndex() - 1; j++) {
            int x = Math.max(j, minInterval);
            if (step2 % (1L << x) == 0) {

                Vector[] config1 = configStorage.getSavedConfig(j);
                Vector[] config0 = configStorage.getSavedConfig(j + 1);
                for (int i = 0; i < config0.length; i++) {
                    IAtom iAtom = atoms.get(i);
                    if (type != null && iAtom.getType() != type) continue;
                    Vector ir0 = config0[i];
                    Vector ir1 = config1[i];
                    Vector ir2 = config2[i];
                    dr01.Ev1Mv2(ir1, ir0);
                    dr12.Ev1Mv2(ir2, ir1);
                    double r01sq = dr01.squared();
                    double r12sq = dr12.squared();
                    dr2Sum[j - 1] += r01sq + r12sq;
                    if (correlationType == CorrelationType.TOTAL) {
                        corSum[j - 1] += dr01.dot(dr12);
                    } else if (correlationType == CorrelationType.MAGNITUDE) {
                        dr1Sum[j - 1] += Math.sqrt(r01sq) + Math.sqrt(r12sq);
                        corSum[j - 1] += Math.sqrt(r01sq * r12sq);
                    } else if (correlationType == CorrelationType.MAG_DOT) {
                        double r01 = Math.sqrt(r01sq), r12 = Math.sqrt(r12sq);
                        dr3Sum[j - 1] += r01 * r01sq + r12 * r12sq;
                        corSum[j - 1] += (r01 + r12) * dr01.dot(dr12);
                    }
                    nSamples[j - 1]++;
                }
            }
        }
    }
}
