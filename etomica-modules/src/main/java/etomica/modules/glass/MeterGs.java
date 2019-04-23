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
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;

/**
 * Meter for tabulation of the atomic Gs.
 *
 * @author Sabry Moustafa
 */
public class MeterGs implements ConfigurationStorage.ConfigurationStorageListener, IDataSource, DataSourceIndependent {

    public enum CorrelationType {TOTAL, PARALLEL, PERPENDICULAR, MAGNITUDE}

    protected final ConfigurationStorage configStorage;
    protected final CorrelationType correlationType;
    protected final DataSourceUniform xDataSource;
    protected final DataTag tag;
    protected final Vector dr, dri, drj, tmp;
    protected long[] gsSum;
    protected DataFunction data;
    protected DataDoubleArray rData;
    protected double xMax;
    protected AtomType type;
    private IDataInfo dataInfo;
    protected int prevSampleIndex;

    public MeterGs(ConfigurationStorage configStorage) {
        this(configStorage, CorrelationType.TOTAL);
    }

    public MeterGs(ConfigurationStorage configStorage, CorrelationType cType) {
        this.configStorage = configStorage;
        this.correlationType = cType;
        Space space = configStorage.getBox().getSpace();
        prevSampleIndex = 0;

        xDataSource = new DataSourceUniform("r", Length.DIMENSION);
        xDataSource.setTypeMax(LimitType.HALF_STEP);
        xDataSource.setTypeMin(LimitType.HALF_STEP);

        rData = (DataDoubleArray) xDataSource.getData();
        data = new DataFunction(new int[]{rData.getLength()});
        gsSum = new long[rData.getLength()];
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

    public void setAtomTypes(AtomType type) {
        this.type = type;
    }

    /**
     * Zero's out the RDF sum tracked by this meter.
     */
    public void reset() {
        rData = (DataDoubleArray) xDataSource.getData();
        xMax = xDataSource.getXMax();
        data = new DataFunction(new int[]{rData.getLength()});
        gsSum = new long[rData.getLength()];
        dataInfo = new DataInfoFunction("c(r)", Null.DIMENSION, this);
        dataInfo.addTag(tag);
        zeroData();
    }

    public void zeroData() {
        for (int i = 0; i < gsSum.length; i++) gsSum[i] = 0;
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

        double gsCount  = 0;
        final double[] y = data.getData();
        double[] r = rData.getData();
        double D = dr.getD();
        for (int i = 0; i < r.length; i++) {
            gsCount += gsSum[i];
        }
        double dx = xMax/rData.getLength();

        for (int i = 0; i < r.length; i++) {
            y[i] = gsSum[i]/(gsCount*dx);
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

        Box box = configStorage.getBox();
        Vector[] config0 = configStorage.getSavedConfig(0);
        Vector[] configPrev = configStorage.getSavedConfig(prevSampleIndex);
        IAtomList atoms = box.getLeafList();

        for (int i = 0; i < config0.length; i++) {
            IAtom iAtom = atoms.get(i);
            if (type != null && (iAtom.getType() != type )) continue; //does not apply if: type is null (total) or iATom is the specified one
            Vector ir0 = config0[i];
            Vector irPrev = configPrev[i];
            dri.Ev1Mv2(ir0, irPrev);
            double r2 = dri.squared();
            if(r2 > xDataSource.getXMax()*xDataSource.getXMax()) continue;
            int index = xDataSource.getIndex(Math.sqrt(r2));  //determine histogram index
            gsSum[index]++;                        //add once for each atom
        }

    }
}