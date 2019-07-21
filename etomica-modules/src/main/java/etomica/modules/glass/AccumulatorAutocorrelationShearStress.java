/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataGroup;
import etomica.units.dimensions.Time;

import java.util.ArrayList;
import java.util.List;

/**
 * Accumulator designed to compute the normalized autocorrelation of the
 * shear stress, assuming zero average.  All data is saved and processed when requested.  The autocorrelation for
 * each component is computed and then averaged together
 */
public class AccumulatorAutocorrelationShearStress extends DataAccumulator implements DataSourceIndependent {

    protected DataFunction data;
    protected DataDoubleArray tData;
    protected DataDoubleArray.DataInfoDoubleArray tDataInfo;
    protected final DataTag tTag;
    protected int nMax;
    protected IDataInfo inputDataInfo;
    protected final List<double[]> savedData;
    protected IData sum, product;
    protected int dim;
    protected final double dt;
    protected int nBlocks;

    public AccumulatorAutocorrelationShearStress(int nMax, double dt) {
        this.dt = dt;
        savedData = new ArrayList<>();
        this.nMax = nMax;
        data = new DataFunction(new int[]{nMax + 1});
        nBlocks = 30;
        tTag = new DataTag();
    }

    public void setNMax(int newNMax) {
        nMax = newNMax;
        data = new DataFunction(new int[]{nMax + 1});
        processDataInfo(inputDataInfo);
    }

    public int getNMax() {
        return nMax;
    }


    @Override
    protected boolean addData(IData data) {
        int n = (dim * dim - dim) / 2 + dim;
        double[] copy = new double[n];
        int i = 0;
        for (int k = 0; k < dim; k++) {
            for (int l = k+1; l < dim; l++) {
                copy[i] = data.getValue(k * dim + l);
                i++;
            }
        }
        savedData.add(copy);
        return true;
    }

    @Override
    public void reset() {
        savedData.clear();
        double[] y = data.getData();
        for (int i = 0; i < nMax; i++) {
            y[i] = Double.NaN;
        }
    }

    @Override
    public IData getData() {
        int n = (dim * dim - dim) / 2; //D=2: n=1; D=3: n=3
        double[] sumij = new double[n];
        int myMax = savedData.size() / 50;
        if (myMax > nMax) myMax = nMax;

        double[] y = data.getData();
        for (int j = 0; j <= myMax; j++) {
            for (int k = 0; k < n; k++) sumij[k] = 0;
            for (int i = 0; i < savedData.size() - j; i++) {
                double[] iSaved = savedData.get(i);
                double[] jSaved = savedData.get(i + j);
                for (int k = 0; k < n; k++) {
                    sumij[k] += iSaved[k] * jSaved[k];
                }
            }
            for (int k = 0; k < n; k++) {
                sumij[k] /= savedData.size() - j;
            }

            double allSum = 0;
            for (int k = 0; k < n; k++) {
                allSum += sumij[k];
            }
            y[j] = allSum / n;
        }
        double y0 = y[0];
        for (int j = 0; j <= myMax; j++) {
            y[j] /= y0;
        }

        return data;
    }

    @Override
    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        this.inputDataInfo = inputDataInfo;
        sum = inputDataInfo.makeData();
        product = inputDataInfo.makeData();
        dim = (int) Math.round(Math.sqrt(inputDataInfo.getLength()));
        tData = new DataDoubleArray(nMax + 1);
        double[] t = tData.getData();
        double[] y = data.getData();
        for (int i = 0; i <= nMax; i++) {
            t[i] = dt * i;
            y[i] = Double.NaN;
        }
        tDataInfo = new DataDoubleArray.DataInfoDoubleArray("lag time", Time.DIMENSION, new int[]{nMax + 1});
        tDataInfo.addTag(tTag);
        dataInfo = new DataFunction.DataInfoFunction(inputDataInfo.getLabel() + " autocorrelation", inputDataInfo.getDimension(), this);
        dataInfo.addTag(tag);
        return dataInfo;
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
        return tag;
    }
}