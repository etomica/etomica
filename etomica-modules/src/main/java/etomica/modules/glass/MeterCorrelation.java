/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.glass;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.DataSourceUniform.LimitType;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;

import java.util.Arrays;

/**
 * Meter for tabulation of the atomic radial distribution function (RDF).  The
 * meter takes data via actionPerformed and returns the average RDF via
 * getData.
 *
 * @author David Kofke
 */
public class MeterCorrelation implements ConfigurationStorage.ConfigurationStorageListener, IDataSource, DataSourceIndependent {

    public enum CorrelationType {TOTAL, PARALLEL, PERPENDICULAR, MAGNITUDE}

    protected final ConfigurationStorage configStorage;
    protected final CorrelationType correlationType;
    protected final DataSourceUniform xDataSource;
    protected final DataTag tag;
    protected final Vector dr, dri, drk, tmp;
    protected double[][] corSum;
    protected long[][] gSum;
    protected DataFunction data;
    protected DataDoubleArray rData;
    protected double xMax;
    protected AtomType type1, type2;
    private IDataInfo dataInfo;
    protected int prevSampleIndex, minPrevSample;
    protected double[] dr2SumA, dr2SumB, dr1SumA, dr1SumB;
    protected long[] dr2CountA, dr2CountB;

    public MeterCorrelation(ConfigurationStorage configStorage) {
        this(configStorage, CorrelationType.TOTAL);
    }

    public MeterCorrelation(ConfigurationStorage configStorage, CorrelationType cType) {
        this.configStorage = configStorage;
        this.correlationType = cType;
        Space space = configStorage.getBox().getSpace();
        prevSampleIndex = 0;

        xDataSource = new DataSourceUniform("r", Length.DIMENSION);
        xDataSource.setTypeMax(LimitType.HALF_STEP);
        xDataSource.setTypeMin(LimitType.HALF_STEP);

        rData = (DataDoubleArray) xDataSource.getData();
        data = new DataFunction(new int[]{rData.getLength()});
        gSum = new long[0][0];
        corSum = new double[0][0];
        dr2SumA = dr2SumB = dr1SumA = dr1SumB = new double[0];
        dr2CountA = dr2CountB = new long[0];
        dataInfo = new DataInfoFunction("c(r)", Null.DIMENSION, this);

        dr = space.makeVector();
        dri = space.makeVector();
        drk = space.makeVector();
        tmp = space.makeVector();
        tag = new DataTag();
        dataInfo.addTag(tag);
    }

    public void setMinPrevSample(int idx) {
        minPrevSample = idx;
    }

    public int getMinPrevSample() {
        return minPrevSample;
    }

    public void setPrevSampleIndex(int idx) {
        prevSampleIndex = idx;
    }

    public int getPrevSampleIndex() {
        return prevSampleIndex;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public void setAtomTypes(AtomType type1, AtomType type2) {
        this.type1 = type1;
        this.type2 = type2;
    }

    protected void reallocate(int n) {
        if (n == gSum.length) return;
        int oldLength = gSum.length;
        gSum = Arrays.copyOf(gSum, n);
        corSum = Arrays.copyOf(corSum, n);
        dr2SumA = Arrays.copyOf(dr2SumA, n);
        dr2SumB = Arrays.copyOf(dr2SumB, n);
        dr1SumA = Arrays.copyOf(dr1SumA, n);
        dr1SumB = Arrays.copyOf(dr1SumB, n);
        dr2CountA = Arrays.copyOf(dr2CountA, n);
        dr2CountB = Arrays.copyOf(dr2CountB, n);
        for (int i = oldLength; i < n; i++) {
            gSum[i] = new long[xDataSource.getNValues()];
            corSum[i] = new double[xDataSource.getNValues()];
        }
    }

    public void zeroData() {
        reallocate(0);
    }

    /**
     * Returns the RDF, averaged over the calls to actionPerformed since the
     * meter was reset or had some parameter changed (xMax or # of bins).
     */
    public IData getData() {
        if (data.getLength() != data.getLength() ||
                xDataSource.getXMax() != xMax) {
            data = new DataFunction(new int[]{xDataSource.getNValues()});
            xMax = xDataSource.getXMax();
            reallocate(0);
            Arrays.fill(data.getData(), Double.NaN);
            return data;
        }
        if (prevSampleIndex >= minPrevSample + gSum.length || prevSampleIndex - minPrevSample < 0) {
            data.E(Double.NaN);
            return data;
        }
        int j = prevSampleIndex - minPrevSample;
        if (dr2SumA[j] == 0) return data;

        final double[] y = data.getData();
        double[] r = rData.getData();
        double norm = Math.sqrt(dr2SumA[j] / dr2CountA[j] * dr2SumB[j] / dr2CountB[j]);
        double avg1A = dr1SumA[j] / dr2CountA[j];
        double avg1B = dr1SumB[j] / dr2CountB[j];
        double D = dr.getD();
        if (correlationType == CorrelationType.PERPENDICULAR) {
            norm *= (D - 1.0) / D;
        } else if (correlationType == CorrelationType.PARALLEL) {
            norm *= 1.0 / D;
        }
        for (int i = 0; i < r.length; i++) {
            if (correlationType == CorrelationType.MAGNITUDE) {
                y[i] = (corSum[j][i] / gSum[j][i] - avg1A * avg1B) /
                        Math.sqrt((dr2SumA[j] / dr2CountA[j] - avg1A * avg1A) *
                                (dr2SumB[j] / dr2CountB[j] - avg1B * avg1B));
            } else {
                y[i] = corSum[j][i] / (gSum[j][i] * norm);
            }
        }
        return data;
    }

    public DataSourceUniform getXDataSource() {
        return xDataSource;
    }

    public DataDoubleArray getIndependentData(int i) {
        return (DataDoubleArray) xDataSource.getData();
    }

    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return (DataInfoDoubleArray) xDataSource.getDataInfo();
    }

    public DataTag getIndependentTag() {
        return xDataSource.getTag();
    }

    public int getIndependentArrayDimension() {
        return 1;
    }

    @Override
    public void newConfigruation() {
        if (xDataSource.getNValues() != data.getLength() ||
                xDataSource.getXMax() != xMax) {
            data = new DataFunction(new int[]{xDataSource.getNValues()});
            xMax = xDataSource.getXMax();
            reallocate(0);
        }
        long step = configStorage.getSavedSteps()[0];
        if (step % (1L << minPrevSample) != 0) return;

        if (configStorage.getLastConfigIndex() < minPrevSample || minPrevSample == 0) return;

        double xMaxSquared = xMax * xMax;
        Box box = configStorage.getBox();
        Boundary boundary = box.getBoundary();
        // iterate over all pairs
        Vector[] config0 = configStorage.getSavedConfig(0);
        IAtomList atoms = box.getLeafList();
        int idxSum = 0, idxProd = 0;
        if (type1 != null) {
            idxSum = type1.getIndex() + type2.getIndex();
            idxProd = type1.getIndex() * type2.getIndex();
        }

        for (int j = 0; j < configStorage.getLastConfigIndex(); j++) {
            int x = Math.max(j, minPrevSample);
            if (step % (1L << x) == 0) {
                if (j >= gSum.length) reallocate(j + 1);

                Vector[] configPrev = configStorage.getSavedConfig(j + 1);
                for (int i = 0; i < config0.length; i++) {
                    IAtom iAtom = atoms.get(i);
                    if (type1 != null && (iAtom.getType() != type1 && iAtom.getType() != type2)) continue;
                    int iTypeIdx = iAtom.getType().getIndex();
                    Vector ir0 = config0[i];
                    Vector irPrev = configPrev[i];
                    dri.Ev1Mv2(ir0, irPrev);
                    double dri2 = dri.squared();
                    if (type1 == null || iAtom.getType() == type1) {
                        dr2SumA[j] += dri2;
                        dr2CountA[j]++;
                        if (correlationType == CorrelationType.MAGNITUDE) dr1SumA[j] += Math.sqrt(dri2);
                    }
                    if (type2 == null || iAtom.getType() == type2) {
                        dr2SumB[j] += dri2;
                        dr2CountB[j]++;
                        if (correlationType == CorrelationType.MAGNITUDE) dr1SumB[j] += Math.sqrt(dri2);
                    }
                    for (int k = i + 1; k < config0.length; k++) {
                        IAtom kAtom = atoms.get(k);
                        int kTypeIdx = kAtom.getType().getIndex();
                        if (type1 != null && (iTypeIdx * kTypeIdx != idxProd || iTypeIdx + kTypeIdx != idxSum))
                            continue;
                        Vector kr0 = config0[k];

                        dr.Ev1Mv2(ir0, kr0);
                        boundary.nearestImage(dr);
                        double r2 = dr.squared();       //compute pair separation
                        if (r2 < xMaxSquared) {
                            Vector krPrev = configPrev[k];
                            drk.Ev1Mv2(kr0, krPrev);
                            int index = xDataSource.getIndex(Math.sqrt(r2));  //determine histogram index
                            gSum[j][index]++;                        //add once for each atom
                            if (correlationType == CorrelationType.TOTAL) {
                                corSum[j][index] += dri.dot(drk);
                            } else if (correlationType == CorrelationType.MAGNITUDE) {
                                corSum[j][index] += Math.sqrt(dri2 * drk.squared());
                            } else {
                                tmp.Ea1Tv1(drk.dot(dr) / r2, dr);
                                if (correlationType == CorrelationType.PERPENDICULAR) {
                                    drk.ME(tmp);
                                } else {
                                    drk.E(tmp);
                                }
                                tmp.Ea1Tv1(dri.dot(dr) / r2, dr);
                                if (correlationType == CorrelationType.PERPENDICULAR) {
                                    tmp.Ev1Mv2(dri, tmp);
                                }
                                corSum[j][index] += tmp.dot(drk);
                            }
                        }
                    }
                }
            }
        }
    }
}
