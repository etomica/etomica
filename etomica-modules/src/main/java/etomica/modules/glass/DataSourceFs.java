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
import etomica.util.Statefull;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;
import java.util.Arrays;

/**
 * Computes the excess kurtosis (alpha2) for the distribution of displacements
 */
public class DataSourceFs implements IDataSource, ConfigurationStorage.ConfigurationStorageListener, DataSourceIndependent, Statefull {

    protected final ConfigurationStorage configStorage;
    protected DataDoubleArray tData;
    protected DataDoubleArray.DataInfoDoubleArray tDataInfo;
    protected DataFunction data;
    protected DataFunction.DataInfoFunction dataInfo;
    protected double[] fsSum;
    protected final DataTag tTag, tag;
    protected long[] nSamples;
    protected Vector dr, q;
    protected AtomType type;
    protected int minInterval = 3;

    public DataSourceFs(ConfigurationStorage configStorage) {
        this.configStorage = configStorage;
        Space space = configStorage.getBox().getSpace();
        fsSum = new double[0];
        nSamples = new long[0];
        tag = new DataTag();
        tTag = new DataTag();
        dr = space.makeVector();
        q = space.makeVector();
        q.setX(0, 7.0);
        reallocate(0);
    }

    public void setQ(Vector q) {
        for (int i = 0; i < q.getD(); i++) {
            this.q.setX(i, q.getX(i));
        }
        reset();
    }

    public Vector getQ() {
        return q;
    }

    public void reset() {
        reallocate(0);
    }

    public void reallocate(int n) {
        fsSum = Arrays.copyOf(fsSum, n);
        nSamples = Arrays.copyOf(nSamples, n);
        data = new DataFunction(new int[]{n});
        tData = new DataDoubleArray(new int[]{n});
        tDataInfo = new DataDoubleArray.DataInfoDoubleArray("t", Time.DIMENSION, new int[]{n});
        dataInfo = new DataFunction.DataInfoFunction("Fs(t)", Null.DIMENSION, this);
        dataInfo.addTag(tag);

        double[] t = tData.getData();
        if (t.length > 0) {
            double dt = configStorage.getDeltaT();
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
        if(type != null){
            Box box = configStorage.getBox();
            nAtoms = box.getNMolecules(type.getSpecies());
        }

        for (int i = 0; i < fsSum.length; i++) {
            y[i] = fsSum[i] / (nAtoms * nSamples[i]) ; // Why subtract "-1" ?
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
        long step = configStorage.getSavedSteps()[0];
        Vector[] positions = configStorage.getSavedConfig(0);
        Box box = configStorage.getBox();
        IAtomList atoms = box.getLeafList();
        for (int i = 0; i < configStorage.getLastConfigIndex(); i++) {
            int x = Math.max(i, minInterval);
            if (step % (1L << x) == 0) {
                if (i >= fsSum.length) reallocate(i + 1);
                Vector[] iPositions = configStorage.getSavedConfig(i + 1);
                for (int j = 0; j < positions.length; j++) {
                    IAtom jAtom = atoms.get(j);
                    if (type == null || jAtom.getType() == type) {
                        dr.Ev1Mv2(positions[j], iPositions[j]);
                        fsSum[i] += Math.cos(q.dot(dr));
                    }
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
        fw.write(fsSum.length+"\n");
        for (int i=0; i<fsSum.length; i++) {
            fw.write(fsSum[i]+" "+nSamples[i]+"\n");
        }
    }

    @Override
    public void restoreState(BufferedReader br) throws IOException {
        int n = Integer.parseInt(br.readLine());
        reallocate(n);
        for (int i=0; i<n; i++) {
            String[] bits = br.readLine().split(" ");
            fsSum[i] = Double.parseDouble(bits[0]);
            nSamples[i] = Long.parseLong(bits[1]);
        }
    }
}
