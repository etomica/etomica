/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.normalmode;

import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.molecule.IMolecule;
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Time;
import etomica.util.Statefull;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;
import java.util.Arrays;

/**
 * Computes the autocorrelation of the vector of ring COM from the lattice site
 */
public class DataSourceRAC implements IDataSource, ConfigurationStorage.ConfigurationStorageListener, DataSourceIndependent, Statefull {

    protected final ConfigurationStorage configStorage;
    protected DataDoubleArray tData;
    protected DataDoubleArray.DataInfoDoubleArray tDataInfo;
    protected DataFunction data, errData;
    protected DataFunction.DataInfoFunction dataInfo;
    protected double[] bacSum, bac2Sum;
    protected final DataTag tTag, tag;
    protected long[] nSamples;
    protected int minInterval = 3;
    protected final Vector[] latticePositions;

    public DataSourceRAC(ConfigurationStorage configStorage) {
        this.configStorage = configStorage;
        bacSum = new double[0];
        bac2Sum = new double[0];
        nSamples = new long[0];
        tag = new DataTag();
        tTag = new DataTag();
        reallocate(0);

        Box box = configStorage.getBox();
        latticePositions = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        for (IMolecule m : box.getMoleculeList()) {
            latticePositions[m.getIndex()].E(m.getChildList().get(0).getPosition());
        }

    }

    public void reset() {
        reallocate(0);
    }

    protected void reallocate(int n) {
        bacSum = Arrays.copyOf(bacSum, n);
        bac2Sum = Arrays.copyOf(bac2Sum, n);
        nSamples = Arrays.copyOf(nSamples, n);
        data = new DataFunction(new int[]{n});
        errData = new DataFunction(new int[]{n});
        tData = new DataDoubleArray(new int[]{n});
        tDataInfo = new DataDoubleArray.DataInfoDoubleArray("t", Time.DIMENSION, new int[]{n});
        dataInfo = new DataFunction.DataInfoFunction("RAC(t)", Null.DIMENSION, this);
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
        double[] yErr = errData.getData();
        int nMolecules = configStorage.getBox().getMoleculeList().size();

        for (int i = 0; i < bacSum.length; i++) {
            long M = nMolecules*nSamples[i];
            y[i] = bacSum[i]/M / (bacSum[0]/(nMolecules*nSamples[0]));
            yErr[i] = Math.sqrt((bac2Sum[i]/M - y[i]*y[i]) / (M - 1)) / (bacSum[0]/(nMolecules*nSamples[0]));
        }
        return data;
    }

    public IData getErrorData() {
        return errData;
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
        Vector[] coords = configStorage.getSavedConfig(0);
        Box box = configStorage.getBox();
        int n = box.getMoleculeList().get(0).getChildList().size();
        for (int i = 0; i < configStorage.getLastConfigIndex(); i++) {
            int x = Math.max(i, minInterval);
            if (step % (1L << x) == 0) {
                if (i >= bacSum.length) reallocate(i + 1);
                Vector[] iCoords = configStorage.getSavedConfig(i + 1);
                for (int j = 0; j < box.getMoleculeList().size(); j++) {
                    double dot = computeDR(coords, j, n).dot(computeDR(iCoords, j, n));
                    bacSum[i] += dot;
                    bac2Sum[i] += dot*dot;
                }
                nSamples[i]++;
            }
        }
    }

    protected Vector computeDR(Vector[] rAtom, int jMol, int n) {
        Vector dr = configStorage.getBox().getSpace().makeVector();
        Vector drSum = configStorage.getBox().getSpace().makeVector();
        Vector site = latticePositions[jMol];
        Boundary boundary = configStorage.getBox().getBoundary();
        for (int i = jMol*n; i < (jMol+1)*n; i++) {
            dr.Ev1Mv2(rAtom[i], site);
            boundary.nearestImage(dr);
            drSum.PE(dr);
        }
        drSum.TE(1.0/n);
        return drSum;
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
        fw.write(bacSum.length+"\n");
        for (int i = 0; i< bacSum.length; i++) {
            fw.write(bacSum[i]+" "+ bac2Sum[i]+" "+nSamples[i]+"\n");
        }
    }

    @Override
    public void restoreState(BufferedReader br) throws IOException {
        int n = Integer.parseInt(br.readLine());
        reallocate(n);
        for (int i=0; i<n; i++) {
            String[] bits = br.readLine().split(" ");
            bacSum[i] = Double.parseDouble(bits[0]);
            bac2Sum[i] = Double.parseDouble(bits[1]);
            nSamples[i] = Long.parseLong(bits[2]);
        }
    }
}