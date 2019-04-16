/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass2d;

import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataGroup;
import etomica.units.dimensions.Time;

import java.util.ArrayList;
import java.util.List;

/**
 * Accumulator designed to compute the unnormalized autocorrelation of the
 * pressure tensor components, assuming that off-diagonal components average to
 * 0.  All data is saved and processed when requested.  The autocorrelation for
 * each component is computed and then averaged together
 */
public class AccumulatorAutocorrelationPTensor extends DataAccumulator implements DataSourceIndependent {

    protected DataFunction data;
    protected DataDoubleArray tData;
    protected DataDoubleArray.DataInfoDoubleArray tDataInfo;
    protected final DataTag tTag, errTag;
    protected int nMax;
    protected IDataInfo inputDataInfo;
    protected final List<double[]> savedData;
    protected IData sum, product;
    protected int dim;
    protected final double dt;
    protected int nBlocks;
    protected DataFunction errData;
    protected DataFunction.DataInfoFunction errDataInfo;
    protected DataFork avgErrFork;
    protected DataGroup avgErrData;
    protected DataGroup.DataInfoGroup avgErrDataInfo;

    public AccumulatorAutocorrelationPTensor(int nMax, double dt) {
        this.dt = dt;
        savedData = new ArrayList<>();
        this.nMax = nMax;
        data = new DataFunction(new int[]{nMax + 1});
        errData = new DataFunction(new int[]{nMax + 1});
        nBlocks = 30;
        tTag = new DataTag();
        errTag = new DataTag();
        avgErrFork = new DataFork();
    }

    public DataFork getAvgErrFork() {
        return avgErrFork;
    }

    @Override
    protected boolean addData(IData data) {
        int n = (dim * dim - dim) / 2 + dim;
        double[] copy = new double[n];
        int i = 0;
        for (int k = 0; k < dim; k++) {
            for (int l = k; l < dim; l++) {
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
        double[] yErr = errData.getData();
        for (int i = 0; i < nMax; i++) {
            y[i] = yErr[i] = Double.NaN;
        }
    }

    @Override
    public IData getData() {
        int n = (dim * dim - dim) / 2 + dim;
        double[] sumij = new double[n];
        double[] sumi = new double[n];
        double[] sumj = new double[n];

        int myMax = savedData.size() / 50;
        if (myMax > nMax) myMax = nMax;

        double[] y = data.getData();
        for (int j = 0; j <= myMax; j++) {
            for (int k = 0; k < n; k++) sumij[k] = sumi[k] = sumj[k] = 0;
            for (int i = 0; i < savedData.size() - j; i++) {
                double[] iSaved = savedData.get(i);
                double[] jSaved = savedData.get(i + j);
                for (int k = 0; k < n; k++) {
                    sumi[k] += iSaved[k];
                    sumj[k] += jSaved[k];
                    sumij[k] += iSaved[k] * jSaved[k];
                }
            }
            for (int k = 0; k < n; k++) {
                sumij[k] /= savedData.size() - j;
                sumi[k] /= savedData.size() - j;
                sumj[k] /= savedData.size() - j;
            }
            int i = 0;
            for (int k = 0; k < dim; k++) {
                sumij[i] -= sumi[i] * sumj[i];
                i += dim - k;
            }
//            System.out.println(j+" sum "+sum[0]+" "+sum[1]+" "+sum[2]);
            double allSum = 0;
            for (int k = 0; k < n; k++) {
                allSum += sumij[k];
            }
            y[j] = allSum / n;
//            System.out.println(j+" "+allSum+" "+n+" "+y[j]);
        }

        int samplesPerBlock = (savedData.size() - myMax) / nBlocks;
        if (samplesPerBlock < 10 || avgErrFork.getDataSinks().length == 0) return data;
        double[] sum1 = new double[myMax + 1];
        double[] sum2 = new double[myMax + 1];
        double[] yErr = errData.getData();
        for (int b = 0; b < nBlocks; b++) {
            for (int j = 0; j <= myMax; j++) {
                for (int k = 0; k < n; k++) sumij[k] = sumi[k] = sumj[k] = 0;
                for (int i = b * samplesPerBlock; i < (b + 1) * samplesPerBlock; i++) {
                    double[] iSaved = savedData.get(i);
                    double[] jSaved = savedData.get(i + j);
                    for (int k = 0; k < n; k++) {
                        sumi[k] += iSaved[k];
                        sumj[k] += jSaved[k];
                        sumij[k] += iSaved[k] * jSaved[k];
                    }
                }
                for (int k = 0; k < n; k++) {
                    sumij[k] /= samplesPerBlock;
                    sumi[k] /= samplesPerBlock;
                    sumj[k] /= samplesPerBlock;
                }
                int i = 0;
                for (int k = 0; k < dim; k++) {
                    // subtract avg squared value for diagonal components
                    sumij[i] -= sumj[i] * sumi[i];
                    i += dim - k;
                }
//                if (b==0) System.out.println(j+" block "+b+" sum "+sum[0]+" "+sum[1]+" "+sum[2]);
                double allSum = 0;
                for (int k = 0; k < n; k++) {
                    allSum += sumij[k];
                }
                sum1[j] += allSum / n;
                sum2[j] += (allSum / n) * (allSum / n);
//                if (b==0) System.out.println(j+" block "+b+" "+allSum+" "+n+" "+sum1[j]);
            }
        }
        for (int j = 0; j <= myMax; j++) {
            yErr[j] = Math.sqrt((nBlocks * sum2[j] - sum1[j] * sum1[j]) / (nBlocks * (nBlocks - 1)));
//            y[j] = sum1[j]/nBlocks;
        }
        avgErrFork.putData(avgErrData);
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
        errDataInfo = new DataFunction.DataInfoFunction(inputDataInfo.getLabel() + " autocorrelation error", inputDataInfo.getDimension(), this);
        errDataInfo.addTag(errTag);
        avgErrData = new DataGroup(new IData[]{data, errData});
        avgErrDataInfo = new DataGroup.DataInfoGroup("autocorrelation", inputDataInfo.getDimension(), new IDataInfo[]{dataInfo, errDataInfo});
        avgErrFork.putDataInfo(avgErrDataInfo);
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
