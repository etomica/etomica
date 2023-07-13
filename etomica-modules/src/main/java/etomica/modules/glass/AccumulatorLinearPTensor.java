/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Time;
import etomica.util.Statefull;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;

public class AccumulatorLinearPTensor implements IDataSink, IDataSource, DataSourceIndependent, Statefull {

    protected DataDoubleArray xData;
    protected DataDoubleArray.DataInfoDoubleArray xDataInfo;
    protected DataFunction data, errData;
    protected DataFunction.DataInfoFunction dataInfo, errDataInfo;
    protected final DataTag xTag, tag, errTag;
    protected double[][] blockSumX;
    protected long[] nBlockSamplesX;
    protected long nSample = 0;
    protected boolean enabled;
    protected double[] sumP2, sumP4;
    protected int nP, dim;
    protected double volume, temperature, dt;
    protected List<IDataSourceCorBlock> viscositySinks;

    public AccumulatorLinearPTensor(Box box, double dt, int maxStored) {
        this.dt = dt;
        this.nP = box.getSpace().D() == 2 ? 1 : 3;
        tag = new DataTag();
        errTag = new DataTag();
        xTag = new DataTag();
        nBlockSamplesX = new long[maxStored];
        blockSumX = new double[maxStored][nP];
        sumP2 = new double[maxStored];
        sumP4 = new double[maxStored];
        volume = box.getBoundary().volume();
        temperature = 1;
        dim = box.getSpace().D();
        viscositySinks = new ArrayList<>();
        reset();
    }

    public void addViscositySink(IDataSourceCorBlock sink) {
        viscositySinks.add(sink);
    }

    public void reset() {
        for (int i = 0; i < nBlockSamplesX.length; i++) {
            nBlockSamplesX[i] = 0;
            sumP2[i] = 0;
            sumP4[i] = 0;
            for(int k=0; k<nP; k++){
                blockSumX[i][k] = 0;
            }
        }

        data = new DataFunction(new int[]{sumP2.length});
        errData = new DataFunction(new int[]{sumP2.length});
        xData = new DataDoubleArray(new int[]{sumP2.length});
        double[] x = xData.getData();
        for(int i=0; i<xData.getLength(); i++){
            x[i] = dt*(i+1);
        }
        xDataInfo = new DataDoubleArray.DataInfoDoubleArray("time", Time.DIMENSION, new int[]{sumP2.length});
        xDataInfo.addTag(xTag);
        dataInfo = new DataFunction.DataInfoFunction("viscosity", Null.DIMENSION, this);
        dataInfo.addTag(tag);
        errDataInfo = new DataFunction.DataInfoFunction("viscosity error", Null.DIMENSION, this);
        errDataInfo.addTag(errTag);
    }


    public void setEnabled(boolean isEnabled) {
        enabled = isEnabled;
        reset();
    }

    public boolean getEnabled() {
        return enabled;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    @Override
    public void putData(IData inputData) {
        if (!enabled) return;

        double[] x = new double[nP];

        for (int i=0, k = 0; k < dim; k++) {
            for (int l = k+1; l < dim; l++) {
                x[i] = inputData.getValue(k * dim + l);
                i++;
            }
        }

        nSample += 1;
        for (int i = 0; i < blockSumX.length; i++) {
            for(int k=0; k<nP; k++){
                blockSumX[i][k] += x[k];
            }
            if (nSample % (i+1) == 0) {
                double sumk = 0;
                for(int k=0; k<nP; k++){
                    double p2 = blockSumX[i][k]*blockSumX[i][k];
                    sumP2[i] += p2;
                    sumP4[i] += p2*p2;
                    blockSumX[i][k] = 0;
                    sumk += p2;
                }
                // sumP2 from each dimension is independent, so we can you the average
                // across dimensions.  using p2 instead of sumk above does not change the
                // resulting average or uncertainty, but does make the uncertainty more
                // precise
                sumk /= nP;
                for (IDataSourceCorBlock s : viscositySinks) {
                    s.putBlock(i, nSample, sumk);
                }

                nBlockSamplesX[i]+=nP;
            }
        }
    }


    public IData getData() {
        double[] y = data.getData();
        double[] yErr = errData.getData();

        double fac = dt*volume/(2*temperature);

        for(int i=0; i<sumP2.length; i++){
            long n = (i+1);
            long M = nBlockSamplesX[i];
            y[i] = fac/M/n*sumP2[i];
            if(M != 0){
                yErr[i] = fac/M/n*Math.sqrt((M*sumP4[i] - sumP2[i]*sumP2[i])/(M-1));
            }
        }
        return data;
    }


    @Override
    public void putDataInfo(IDataInfo inputDataInfo) {
    }

    @Override
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    @Override
    public DataDoubleArray getIndependentData(int i) {
        return xData;
    }

    @Override
    public DataDoubleArray.DataInfoDoubleArray getIndependentDataInfo(int i) {
        return xDataInfo;
    }

    @Override
    public int getIndependentArrayDimension() {
        return 1;
    }

    @Override
    public DataTag getIndependentTag() {
        return xTag;
    }


    public TemperatureSink makeTemperatureSink() {
        return new TemperatureSink();
    }

    @Override
    public void saveState(Writer fw) throws IOException {
        fw.write(getClass().getName()+"\n");
        fw.write(""+nSample+"\n");
        if (nSample == 0) {
            return;
        }
        for (int i=0; i<nBlockSamplesX.length && nBlockSamplesX[i] != 0; i++) {
            fw.write(nBlockSamplesX[i]+" "+sumP2[i]+" "+sumP4[i]);
            for (int j=0; j<nP; j++) {
                fw.write(" "+blockSumX[i][j]);
            }
            fw.write("\n");
        }
        fw.write("\n");
    }

    @Override
    public void restoreState(BufferedReader br) throws IOException {
        if (!br.readLine().equals(getClass().getName())) throw new RuntimeException("oops");
        nSample = Long.parseLong(br.readLine());
        if (nSample == 0) return;
        for (int i=0; i<nBlockSamplesX.length; i++) {
            String s = br.readLine();
            if (s.length() < 2) return;
            String[] bits = s.split(" ");
            nBlockSamplesX[i] = Long.parseLong(bits[0]);
            sumP2[i] = Double.parseDouble(bits[1]);
            sumP4[i] = Double.parseDouble(bits[2]);
            for (int j=0; j<nP; j++) {
                blockSumX[i][j] = Double.parseDouble(bits[3+j]);
            }
        }
    }

    public class TemperatureSink implements IDataSink {
        public void putDataInfo(IDataInfo inputDataInfo) {
        }

        public void putData(IData data) {
            temperature = data.getValue(0);
        }
    }

    public interface ViscositySink {
        void putViscosity(int log2interval, long step, double viscosity);
    }
}