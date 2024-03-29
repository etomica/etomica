/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.normalmode;

import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.space.Vector;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Time;
import etomica.util.Statefull;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;
import java.util.Arrays;

/**
 * Computes the autocorrelation of the bond vector
 */
public class DataSourceBAC implements IDataSource, ConfigurationStorage.ConfigurationStorageListener, DataSourceIndependent, Statefull {

    protected final ConfigurationStorage configStorage;
    protected DataDoubleArray tData;
    protected DataDoubleArray.DataInfoDoubleArray tDataInfo;
    protected DataFunction data, errData;
    protected DataFunction.DataInfoFunction dataInfo;
    protected double[] bacSum, bac2Sum;
    protected final DataTag tTag, tag;
    protected long[] nSamples;
    protected int minInterval = 3;
    protected double sum0;

    public DataSourceBAC(ConfigurationStorage configStorage) {
        this.configStorage = configStorage;
        bacSum = new double[0];
        bac2Sum = new double[0];
        nSamples = new long[0];
        tag = new DataTag();
        tTag = new DataTag();
        reallocate(0);
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
        dataInfo = new DataFunction.DataInfoFunction("BAC(t)", Null.DIMENSION, this);
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
        int nAtoms = configStorage.getSavedConfig(0).length;
        double norm = sum0/(nAtoms*nSamples[0]);

        for (int i = 0; i < bacSum.length; i++) {
            long M = nAtoms*nSamples[i];
            y[i] = bacSum[i]/M;
            yErr[i] = Math.sqrt((bac2Sum[i]/M - y[i]*y[i]) / (M - 1)) / norm;
            y[i] /= norm;
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
        Vector dr = box.getSpace().makeVector();
        Vector idr = box.getSpace().makeVector();
        for (int i = 0; i < configStorage.getLastConfigIndex(); i++) {
            int x = Math.max(i, minInterval);
            if (step % (1L << x) == 0) {
                if (i >= bacSum.length) reallocate(i + 1);
                Vector[] iCoords = configStorage.getSavedConfig(i + 1);
                for (int j = 0; j < coords.length; j++) {
                    int jp1 = j+1;
                    if (jp1%n == 0) jp1-=n;
                    dr.Ev1Mv2(coords[j], coords[jp1]);
                    box.getBoundary().nearestImage(dr);
                    idr.Ev1Mv2(iCoords[j], iCoords[jp1]);
                    box.getBoundary().nearestImage(idr);
                    double dot = dr.dot(idr);
                    bacSum[i] += dot;
                    bac2Sum[i] += dot*dot;
                    if (i == 0) sum0 += dr.squared();
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