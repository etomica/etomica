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
import etomica.util.Statefull;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;
import java.util.Arrays;

/**
 * Meter for tabulation of the atomic Gs.
 *
 * @author Sabry Moustafa
 */
public class MeterGs implements ConfigurationStorage.ConfigurationStorageListener, IDataSource, DataSourceIndependent, Statefull {

    protected final ConfigurationStorage configStorage;
    protected final DataSourceUniform xDataSource;
    protected final DataTag tag;
    protected final Vector dr, dri, drj, tmp;
    protected int minConfigIndex, configIndex;
    protected long[][] gsSum;
    protected DataFunction data;
    protected double xMax;
    protected AtomType type;
    protected IDataInfo dataInfo;
    protected long[] numSamples;

    public MeterGs(ConfigurationStorage configStorage) {
        this.configStorage = configStorage;

        xDataSource = new DataSourceUniform("r", Length.DIMENSION);
        xDataSource.setTypeMax(LimitType.HALF_STEP);
        xDataSource.setTypeMin(LimitType.HALF_STEP);

        gsSum = new long[0][0];
        data = new DataFunction(new int[]{xDataSource.getNValues()});
        tag = new DataTag();
        dataInfo = new DataInfoFunction("c(r)", Null.DIMENSION, this);
        dataInfo.addTag(tag);

        numSamples = new long[0];

        Space space = configStorage.getBox().getSpace();
        dr = space.makeVector();
        dri = space.makeVector();
        drj = space.makeVector();
        tmp = space.makeVector();
    }

    protected void reallocate(int n) {
        if (n == gsSum.length) return;
        int oldLength = gsSum.length;
        gsSum = Arrays.copyOf(gsSum, n);
        numSamples = Arrays.copyOf(numSamples, n);
        for (int i = oldLength; i < n; i++) {
            gsSum[i] = new long[xDataSource.getNValues()];
        }
    }

    public void reset() {
        reallocate(0);
    }

    public void setMinConfigIndex(int idx) {
        minConfigIndex = idx;
    }

    public int getMinConfigIndex() {
        return minConfigIndex;
    }

    public void setConfigIndex(int idx) {
        configIndex = idx;
    }

    public int getConfigIndex() {
        return configIndex;
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
        if (configIndex >= gsSum.length) {
            data.E(Double.NaN);
            return data;
        }

        double gsCount = numSamples[configIndex];
        final double[] y = data.getData();
        double[] r = ((DataDoubleArray) xDataSource.getData()).getData();
        double dx = xMax / xDataSource.getNValues();

        for (int i = 0; i < r.length; i++) {
            y[i] = gsSum[configIndex][i] / (gsCount * dx);
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
        if (step % (1L << minConfigIndex) != 0) return;

        Box box = configStorage.getBox();
        Vector[] config0 = configStorage.getSavedConfig(0);
        IAtomList atoms = box.getLeafList();
        double xMax2 = xMax * xMax;

        for (int j = 0; j < configStorage.getLastConfigIndex(); j++) {
            int x = Math.max(j, minConfigIndex);
            if (step % (1L << x) == 0) {
                if (j >= gsSum.length) reallocate(j + 1);
                Vector[] configPrev = configStorage.getSavedConfig(j + 1);

                for (int i = 0; i < config0.length; i++) {
                    IAtom iAtom = atoms.get(i);
                    if (type != null && (iAtom.getType() != type))
                        continue; //does not apply if: type is null (total) or iATom is the specified one
                    numSamples[j]++;
                    Vector ir0 = config0[i];
                    Vector irPrev = configPrev[i];
                    dri.Ev1Mv2(ir0, irPrev);
                    double r2 = dri.squared();
                    if (r2 > xMax2) continue;
                    int index = xDataSource.getIndex(Math.sqrt(r2));  //determine histogram index
                    gsSum[j][index]++;                        //add once for each atom
                }
            }
        }
    }

    @Override
    public void saveState(Writer fw) throws IOException {
        fw.write(gsSum.length+"\n");
        for (int i=0; i<gsSum.length; i++) {
            fw.write(""+numSamples[i]);
            for (int j=0; j<gsSum[i].length; j++) {
                fw.write(" "+gsSum[i][j]);
            }
            fw.write("\n");
        }
    }

    @Override
    public void restoreState(BufferedReader br) throws IOException {
        int n = Integer.parseInt(br.readLine());
        reallocate(n);
        for (int i=0; i<n; i++) {
            String[] bits = br.readLine().split(" ");
            numSamples[i] = Long.parseLong(bits[0]);
            for (int j=0; j<gsSum[i].length; j++) {
                gsSum[i][j] = Long.parseLong(bits[1+j]);
            }
        }
    }
}
