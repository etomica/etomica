/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
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
public class DataSourceVAC implements IDataSource, ConfigurationStorage.ConfigurationStorageListener, DataSourceIndependent {

    protected final ConfigurationStorage configStorage;
    protected DataDoubleArray tData;
    protected DataDoubleArray.DataInfoDoubleArray tDataInfo;
    protected DataFunction data, errData;
    protected DataFunction.DataInfoFunction dataInfo;
    protected double[] vacSum, vac2Sum;
    protected final DataTag tTag, tag;
    protected long[] nSamples;
    protected AtomType type;
    protected Space space;

    public DataSourceVAC(ConfigurationStorage configStorage) {
        this.configStorage = configStorage;
        space = configStorage.getBox().getSpace();
        vacSum = new double[0];
        vac2Sum = new double[0];
        nSamples = new long[0];
        tag = new DataTag();
        tTag = new DataTag();
        reset();
    }


    public void reset() {
        int n = configStorage.getLastConfigIndex();
        if (n  == vacSum.length && data != null) return;
        if (n < 1) n = 0;
        vacSum = Arrays.copyOf(vacSum, n);
        vac2Sum = Arrays.copyOf(vac2Sum, n);
        nSamples = Arrays.copyOf(nSamples, n);
        data = new DataFunction(new int[]{n});
        errData = new DataFunction(new int[]{n});
        tData = new DataDoubleArray(new int[]{n});
        tDataInfo = new DataDoubleArray.DataInfoDoubleArray("t", Time.DIMENSION, new int[]{n});
        dataInfo = new DataFunction.DataInfoFunction("VAC(t)", Null.DIMENSION, this);
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
        double[] yErr = errData.getData();
        int nAtoms = configStorage.getSavedConfig(0).length;
        if(type != null){
            Box box = configStorage.getBox();
            nAtoms = box.getNMolecules(type.getSpecies());
        }

        for (int i = 0; i < vacSum.length; i++) {
            long M = nAtoms*nSamples[i];
            y[i] = vacSum[i]/M;
            yErr[i] = Math.sqrt((vac2Sum[i]/M - y[i]*y[i]) / (M - 1));
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

    public void setAtomType(AtomType type) {
        this.type = type;
    }

    @Override
    public void newConfigruation() {
        reset(); // reallocates if needed
        long step = configStorage.getSavedSteps()[0];
        Vector[] velocities = configStorage.getSavedVel(0);
        Box box = configStorage.getBox();
        IAtomList atoms = box.getLeafList();
        for (int i = 1; i < vacSum.length; i++) {
            if (step % (1L << (i - 1)) == 0) {
                Vector[] iVelocities = configStorage.getSavedVel(i);
                for (int j = 0; j < velocities.length; j++) {
                    IAtom jAtom = atoms.get(j);
                    if(type == null || jAtom.getType() == type){
                        double vaci = velocities[j].dot(iVelocities[j]);
                        vacSum[i-1] += vaci;
                        vac2Sum[i-1] += vaci*vaci;
                    }
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