/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.molecule.IMoleculeList;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Time;
import etomica.util.Statefull;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;
import java.util.Arrays;

/**
 * Computes the excess kurtosis (alpha2) for the distribution of displacements
 */
public class DataSourceAlpha2 implements IDataSource, ConfigurationStorage.ConfigurationStorageListener, DataSourceIndependent, Statefull {

    protected final ConfigurationStorage configStorage;
    protected DataDoubleArray tData;
    protected DataDoubleArray.DataInfoDoubleArray tDataInfo;
    protected DataFunction data;
    protected DataFunction.DataInfoFunction dataInfo;
    protected double[] msdSum, m4dSum;
    protected final DataTag tTag, tag;
    protected long[] nSamples;
    protected int minInterval = 4;
    protected final ISpecies species;

    public DataSourceAlpha2(ConfigurationStorage configStorage) {
        this(configStorage, null);
    }

    public DataSourceAlpha2(ConfigurationStorage configStorage, ISpecies species) {
        this.species = species;
        this.configStorage = configStorage;
        msdSum = new double[0];
        m4dSum = new double[0];
        nSamples = new long[0];
        tag = new DataTag();
        tTag = new DataTag();
        reallocate(0);
    }

    public void reallocate(int n) {
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
            double dt = configStorage.getDeltaT();
            for (int i = 0; i < t.length; i++) {
                t[i] = dt * (1L << i);
            }
        }
    }

    public void reset() {
        reallocate(0);
    }

    @Override
    public IData getData() {
        if (configStorage.getLastConfigIndex() < 1) return data;
        double[] y = data.getData();
        int nAtoms = configStorage.getSavedConfig(0).length;
        // (3/5) for 3D; (1/2) for 2D
        double fac = configStorage.getBox().getSpace().D() == 2 ? 0.5 : 0.6;

        if (species != null) {
            nAtoms = configStorage.getBox().getNMolecules(species);
        }

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
        long step = configStorage.getSavedSteps()[0];
        Vector[] positions = configStorage.getSavedConfig(0);
        IMoleculeList molecules = configStorage.getBox().getMoleculeList();
        for (int i = 0; i < configStorage.getLastConfigIndex(); i++) {
            int x = Math.max(i, minInterval);
            if (step % (1L << x) == 0) {
                if (i >= msdSum.length) reallocate(i + 1);
                Vector[] iPositions = configStorage.getSavedConfig(i + 1);
                for (int j = 0; j < positions.length; j++) {
                    if (species != null && molecules.get(j).getType() != species) continue;
                    double d2 = positions[j].Mv1Squared(iPositions[j]);
                    msdSum[i] += d2;
                    m4dSum[i] += d2 * d2;
                }
                nSamples[i]++;
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

    @Override
    public void saveState(Writer fw) throws IOException {
        fw.write(nSamples.length+"\n");
        for (int i=0; i<nSamples.length; i++) {
            fw.write(msdSum[i]+" "+m4dSum[i]+" "+nSamples[i]+"\n");
        }
    }

    @Override
    public void restoreState(BufferedReader br) throws IOException {
        int n = Integer.parseInt(br.readLine());
        reallocate(n);
        for (int i=0; i<n; i++) {
            String[] bits = br.readLine().split(" ");
            msdSum[i] = Double.parseDouble(bits[0]);
            m4dSum[i] = Double.parseDouble(bits[1]);
            nSamples[i] = Long.parseLong(bits[2]);
        }
    }
}
