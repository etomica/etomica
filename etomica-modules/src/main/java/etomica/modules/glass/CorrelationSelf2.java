/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.glass;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.data.DataSourceIndependentSimple;
import etomica.data.DataTag;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;

import java.util.Arrays;

public class CorrelationSelf2 implements ConfigurationStorage.ConfigurationStorageListener {

    public enum CorrelationType {TOTAL, MAGNITUDE}

    protected CorrelationType correlationType;
    protected final ConfigurationStorage configStorage;
    protected final DataTag[] tag, rTag;
    protected final Vector dr01, dr12;
    protected final double[][] corSum, dr2Sum;
    protected final long[][] nSamples;
    protected final DataFunction[] data;
    protected final DataDoubleArray[] rData;
    protected AtomType type;
    protected final DataFunction.DataInfoFunction[] dataInfo;
    protected double[] tdr;
    protected int[] iminDr;
    protected final double dr0;
    protected int minInterval = 3;

    public CorrelationSelf2(ConfigurationStorage configStorage, CorrelationType cType, double dr0, int ndr) {
        this.configStorage = configStorage;
        this.correlationType = cType;
        this.dr0 = dr0;
        int ndt = 60;
        data = new DataFunction[ndt];
        rData = new DataDoubleArray[ndt];
        dataInfo = new DataFunction.DataInfoFunction[ndt];
        rTag = new DataTag[ndt];
        tag = new DataTag[ndt];
        for (int i = 0; i < ndt; i++) {
            data[i] = new DataFunction(new int[]{ndr});
            rData[i] = new DataDoubleArray(new int[]{ndr});
            DataDoubleArray.DataInfoDoubleArray rDataInfo = new DataDoubleArray.DataInfoDoubleArray("r", Length.DIMENSION, new int[]{ndr});
            dataInfo[i] = new DataFunction.DataInfoFunction("cor", Null.DIMENSION, new DataSourceIndependentSimple(rData[i].getData(), rDataInfo));
            tag[i] = new DataTag();
            dataInfo[i].addTag(tag[i]);
            rTag[i] = new DataTag();
        }
        Space space = configStorage.getBox().getSpace();
        dr2Sum = new double[ndr][0];
        corSum = new double[ndr][0];
        nSamples = new long[ndr][0];

        dr01 = space.makeVector();
        dr12 = space.makeVector();
        iminDr = new int[0];
        tdr = new double[0];
        reallocate();
    }

    public int getNumDr() {
        return corSum.length;
    }

    public int getNumDt() {
        return data.length;
    }

    public void reallocate() {
        int n = configStorage.getLastConfigIndex();
        if (n + 1 == corSum.length && data != null) return;
        if (n < 1) n = 0;
        else n--;
        for (int i = 0; i < corSum.length; i++) {
            corSum[i] = Arrays.copyOf(corSum[i], n);
            nSamples[i] = Arrays.copyOf(nSamples[i], n);
            dr2Sum[i] = Arrays.copyOf(dr2Sum[i], n);
        }
        iminDr = Arrays.copyOf(iminDr, n);
        tdr = Arrays.copyOf(tdr, n);
    }

    public void setAtomType(AtomType type) {
        this.type = type;
    }

    @Override
    public void newConfigruation() {

        reallocate(); // reallocates if needed
        long step2 = configStorage.getSavedSteps()[0];
        Vector[] config2 = configStorage.getSavedConfig(0);
        IAtomList atoms = configStorage.getBox().getLeafList();

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
                    double r01 = Math.sqrt(r01sq);
                    double r12 = Math.sqrt(r12sq);
                    if (tdr[j - 1] == 0) {
                        // first sample.  use dr0 and start range as far to the left (closest to 0) as possible
                        tdr[j - 1] = dr0;
                        iminDr[j - 1] = 0;
                        if (r01 / tdr[j - 1] > corSum.length) {
                            iminDr[j - 1] = (int) Math.ceil(r01 / tdr[j - 1]) - corSum.length;
                        }
                        if (r12 / tdr[j - 1] - iminDr[j - 1] > corSum.length) {
                            iminDr[j - 1] = (int) Math.ceil(r12 / tdr[j - 1]) - corSum.length;
                        }
                    }

                    double dot = dr01.dot(dr12);
                    ensureRange(j - 1, r01);
                    int ir01 = (int) (r01 / tdr[j - 1]) - iminDr[j - 1];
                    corSum[ir01][j - 1] += dot;
                    nSamples[ir01][j - 1]++;
                    dr2Sum[ir01][j - 1] += r01sq + r12sq;

                    ensureRange(j - 1, r12);
                    int ir12 = (int) (r12 / tdr[j - 1]) - iminDr[j - 1];
                    corSum[ir12][j - 1] += dot;
                    nSamples[ir12][j - 1]++;
                    dr2Sum[ir12][j - 1] += r01sq + r12sq;
                }
            }
        }
    }

    /**
     * Collapses data by 2.  The new data range ends at the same position as the old.
     * The new min might be negative.
     */
    protected void collapse(int idt) {
        int n = corSum.length;
        int shift = iminDr[idt] % 2;
        if (shift == 0) {
            // for shift=1, the last bin is still just the last bin
            corSum[n - 1][idt] += corSum[n - 2][idt];
            dr2Sum[n - 1][idt] += dr2Sum[n - 2][idt];
            nSamples[n - 1][idt] += nSamples[n - 2][idt];
        }
        for (int i = 1; i < n; i++) {
            int foo = n - 1 - 2 * i + shift;
            if (foo < 0) {
                // this i covers range not included before
                corSum[n - 1 - i][idt] = 0;
                dr2Sum[n - 1 - i][idt] = 0;
                nSamples[n - 1 - i][idt] = 0;
            } else if (foo == 0) {
                // half of i's new range was included before
                corSum[n - 1 - i][idt] = corSum[n - 1 - 2 * i + shift][idt];
                dr2Sum[n - 1 - i][idt] = dr2Sum[n - 1 - 2 * i + shift][idt];
                nSamples[n - 1 - i][idt] = nSamples[n - 1 - 2 * i + shift][idt];
            } else {
                corSum[n - 1 - i][idt] = corSum[n - 1 - 2 * i + shift][idt] + corSum[n - 1 - 2 * i + shift - 1][idt];
                dr2Sum[n - 1 - i][idt] = dr2Sum[n - 1 - 2 * i + shift][idt] + dr2Sum[n - 1 - 2 * i + shift - 1][idt];
                nSamples[n - 1 - i][idt] = nSamples[n - 1 - 2 * i + shift][idt] + nSamples[n - 1 - 2 * i + shift - 1][idt];
            }
        }
        iminDr[idt] = (iminDr[idt] + n + shift) / 2 - n;
        tdr[idt] *= 2;
        // negative iminDr here isn't actually a problem
    }

    /**
     * Ensures that the range for the given time interval includes the given r
     */
    protected void ensureRange(int idt, double r) {
        boolean changed = false;
        int ir = (int) (r / tdr[idt]) - iminDr[idt];
        while (ir < 0) {
            // our new r is to the left of our range and we are already shifted to the
            // left as much as possible, so we have to collapse.
            collapse(idt);
            ir = (int) (r / tdr[idt]) - iminDr[idt];
            changed = true;
        }
        while (ir >= corSum.length) {
            int min;
            for (min = 0; min < nSamples.length && nSamples[min][idt] == 0; min++) {
            }
            if (ir - min < corSum.length) {
                // shift alone is sufficient
                int shift = 1 + ir - corSum.length;
                for (int k = 0; k < corSum.length - shift; k++) {
                    corSum[k][idt] = corSum[k + shift][idt];
                    dr2Sum[k][idt] = dr2Sum[k + shift][idt];
                    nSamples[k][idt] = nSamples[k + shift][idt];
                }
                for (int k = corSum.length - shift; k < corSum.length; k++) {
                    dr2Sum[k][idt] = corSum[k][idt] = nSamples[k][idt] = 0;
                }
                iminDr[idt] += shift;
            } else {
                // need to collapse.  we produce new bins on the left, which
                // won't immediately help, but we might be able to shift next time
                collapse(idt);
            }
            ir = (int) (r / tdr[idt]) - iminDr[idt];
            changed = true;
        }
        if (!changed) return;

        // now update rData
        double[] x = rData[idt].getData();
        for (int i = 0; i < x.length; i++) {
            x[i] = (iminDr[idt] + i + 0.5) * tdr[idt];
        }
    }

    public MeterCorrelationSelf2 makeMeter(int idr) {
        return new MeterCorrelationSelf2(idr);
    }

    public class MeterCorrelationSelf2 implements IDataSource {

        protected final int mydt;

        public MeterCorrelationSelf2(int idt) {
            mydt = idt;
        }

        @Override
        public DataFunction getData() {
            DataFunction myData = data[mydt];
            if (configStorage.getLastConfigIndex() < Math.max(2, mydt + 2) || mydt >= corSum[0].length) {
                myData.E(Double.NaN);
                return myData;
            }

            final double[] y = myData.getData();
            for (int i = 0; i < y.length; i++) {
                y[i] = 2 * corSum[i][mydt] / dr2Sum[i][mydt];
            }
            return myData;
        }

        @Override
        public DataTag getTag() {
            return tag[mydt];
        }

        @Override
        public DataFunction.DataInfoFunction getDataInfo() {
            return dataInfo[mydt];
        }
    }
}
