/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.glass2d;

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
    protected final Vector dr, dri, drj, tmp;
    protected double[] corSum;
    protected long[] gSum;
    protected DataFunction data;
    protected DataDoubleArray rData;
    protected double xMax;
    protected AtomType type1, type2;
    private IDataInfo dataInfo;
    protected int prevSampleIndex;
    protected double dr2SumA, dr2SumB, dr1SumA, dr1SumB;
    protected long dr2CountA, dr2CountB;

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
        gSum = new long[rData.getLength()];
        corSum = new double[rData.getLength()];
        dataInfo = new DataInfoFunction("c(r)", Null.DIMENSION, this);

        dr = space.makeVector();
        dri = space.makeVector();
        drj = space.makeVector();
        tmp = space.makeVector();
        tag = new DataTag();
        dataInfo.addTag(tag);
    }

    public void setPrevSampleIndex(int idx) {
        prevSampleIndex = idx;
        reset();
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

//    public void setPrevConfigIndex(int newPrevConfigIndex) {
//        prevConfigIndex = newPrevConfigIndex;
//        reset();
//    }
//
//    public int getPrevConfigIndex() {
//        return prevConfigIndex;
//    }

    /**
     * Zero's out the RDF sum tracked by this meter.
     */
    public void reset() {
        rData = (DataDoubleArray) xDataSource.getData();
        xMax = xDataSource.getXMax();
        data = new DataFunction(new int[]{rData.getLength()});
        gSum = new long[rData.getLength()];
        corSum = new double[rData.getLength()];
        dataInfo = new DataInfoFunction("c(r)", Null.DIMENSION, this);
        dataInfo.addTag(tag);
        zeroData();
    }

    public void zeroData() {
        for (int i = 0; i < gSum.length; i++) corSum[i] = gSum[i] = 0;
        dr1SumA = dr1SumB = dr2SumA = dr2SumB = 0;
        dr2CountA = dr2CountB = 0;
    }

    /**
     * Returns the RDF, averaged over the calls to actionPerformed since the
     * meter was reset or had some parameter changed (xMax or # of bins).
     */
    public IData getData() {
        if (rData != xDataSource.getData() ||
                data.getLength() != rData.getLength() ||
                xDataSource.getXMax() != xMax) {
            reset();
            //that zeroed everything.  just return the zeros.
            return data;
        }
        if (dr2SumA == 0) return data;

        final double[] y = data.getData();
        double[] r = rData.getData();
        double norm = Math.sqrt(dr2SumA / dr2CountA * dr2SumB / dr2CountB);
        double avg1A = dr1SumA / dr2CountA;
        double avg1B = dr1SumB / dr2CountB;
        double D = dr.getD();
        if (correlationType == CorrelationType.PERPENDICULAR) {
            norm *= (D - 1.0) / D;
        } else if (correlationType == CorrelationType.PARALLEL) {
            norm *= 1.0 / D;
        }
        for (int i = 0; i < r.length; i++) {
            if (correlationType == CorrelationType.MAGNITUDE) {
                y[i] = (corSum[i] / gSum[i] - avg1A * avg1B) /
                        Math.sqrt((dr2SumA / dr2CountA - avg1A * avg1A) *
                                (dr2SumB / dr2CountB - avg1B * avg1B));
            } else {
                y[i] = corSum[i] / (gSum[i] * norm);
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
        if (rData != xDataSource.getData() ||
                data.getLength() != rData.getLength() ||
                xDataSource.getXMax() != xMax) {
            reset();
        }
        if (configStorage.getLastConfigIndex() < prevSampleIndex || prevSampleIndex == 0) return;

        double xMaxSquared = xMax * xMax;
        Box box = configStorage.getBox();
        Boundary boundary = box.getBoundary();
        // iterate over all pairs
        Vector[] config0 = configStorage.getSavedConfig(0);
        Vector[] configPrev = configStorage.getSavedConfig(prevSampleIndex);
        IAtomList atoms = box.getLeafList();
        int idxSum = 0, idxProd = 0;
        if (type1 != null) {
            idxSum = type1.getIndex() + type2.getIndex();
            idxProd = type1.getIndex() * type2.getIndex();
        }
        for (int i = 0; i < config0.length; i++) {
            IAtom iAtom = atoms.get(i);
            if (type1 != null && (iAtom.getType() != type1 && iAtom.getType() != type2)) continue;
            int iTypeIdx = iAtom.getType().getIndex();
            Vector ir0 = config0[i];
            Vector irPrev = configPrev[i];
            dri.Ev1Mv2(ir0, irPrev);
            double dri2 = dri.squared();
            if (type1 == null || iAtom.getType() == type1) {
                dr2SumA += dri2;
                dr2CountA++;
                if (correlationType == CorrelationType.MAGNITUDE) dr1SumA += Math.sqrt(dri2);
            }
            if (type2 == null || iAtom.getType() == type2) {
                dr2SumB += dri2;
                dr2CountB++;
                if (correlationType == CorrelationType.MAGNITUDE) dr1SumB += Math.sqrt(dri2);
            }
            for (int j = i + 1; j < config0.length; j++) {
                IAtom jAtom = atoms.get(j);
                int jTypeIdx = jAtom.getType().getIndex();
                if (type1 != null && (iTypeIdx * jTypeIdx != idxProd || iTypeIdx + jTypeIdx != idxSum)) continue;
                Vector jr0 = config0[j];

                dr.Ev1Mv2(ir0, jr0);
                boundary.nearestImage(dr);
                double r2 = dr.squared();       //compute pair separation
                if (r2 < xMaxSquared) {
                    Vector jrPrev = configPrev[j];
                    drj.Ev1Mv2(jr0, jrPrev);
                    int index = xDataSource.getIndex(Math.sqrt(r2));  //determine histogram index
                    gSum[index]++;                        //add once for each atom
                    if (correlationType == CorrelationType.TOTAL) {
                        corSum[index] += dri.dot(drj);
                    } else if (correlationType == CorrelationType.MAGNITUDE) {
                        corSum[index] += Math.sqrt(dri2 * drj.squared());
                    } else {
                        tmp.Ea1Tv1(drj.dot(dr) / r2, dr);
                        if (correlationType == CorrelationType.PERPENDICULAR) {
                            drj.ME(tmp);
                        } else {
                            drj.E(tmp);
                        }
                        tmp.Ea1Tv1(dri.dot(dr) / r2, dr);
                        if (correlationType == CorrelationType.PERPENDICULAR) {
                            tmp.Ev1Mv2(dri, tmp);
                        }
                        corSum[index] += tmp.dot(drj);
                    }
                }
            }
        }
    }
}
