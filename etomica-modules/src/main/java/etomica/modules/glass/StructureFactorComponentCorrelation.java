package etomica.modules.glass;

import etomica.data.*;
import etomica.data.meter.MeterStructureFactor;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.space.Vector;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Time;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class StructureFactorComponentCorrelation implements DataSourceIndependent, DataSinkBlockAveragerSFac.Sink {

    protected double[][][] lastXY;
    protected double[][] corSum;
    protected double[][] sum2;
    protected DataDoubleArray tData;
    protected DataDoubleArray.DataInfoDoubleArray tDataInfo;
    protected final DataFunction[] data;
    protected final DataFunction.DataInfoFunction[] dataInfo;
    protected final DataTag[] tag;
    protected final DataTag tTag;
    protected final ConfigurationStorage configStorage;
    protected long[][] nSamplesCor, nSamples2;
    protected int minInterval;
    protected boolean[] skipNow;
    protected final int[] wvMap;

    public StructureFactorComponentCorrelation(int nWaveVectors, ConfigurationStorage configStorage) {
        this(makeIdentityMap(nWaveVectors), configStorage);
    }

    public StructureFactorComponentCorrelation(int[] wvMap, ConfigurationStorage configStorage) {
        this.wvMap = wvMap;
        int nWaveVectors = wvMap.length;
        this.configStorage = configStorage;
        data = new DataFunction[nWaveVectors];
        dataInfo = new DataFunction.DataInfoFunction[nWaveVectors];
        tag = new DataTag[nWaveVectors];
        for (int i = 0; i < nWaveVectors; i++) {
            tag[i] = new DataTag();
        }
        tTag = new DataTag();
        lastXY = new double[nWaveVectors][0][0];
        int maxWV = 0;
        for (int i : wvMap) {
            if (i > maxWV) maxWV = i;
        }
        maxWV++;
        corSum = new double[maxWV][0];
        sum2 = new double[maxWV][0];
        nSamples2 = new long[maxWV][0];
        nSamplesCor = new long[maxWV][0];
        skipNow = new boolean[0];
        reallocate(0);
    }

    private static int[] makeIdentityMap(int n) {
        int[] map = new int[n];
        for (int i = 0; i < n; i++) map[i] = i;
        return map;
    }

    public static int[] makeWaveVectorMap(Vector[] waveVectors, int motionDim) {
        int[] map = new int[waveVectors.length];
        List<Double> wv2List = new ArrayList<>();
        List<Double> wv1List = new ArrayList<>();
        for (int i = 0; i < waveVectors.length; i++) {
            map[i] = -1;
            if (motionDim < 0) {
                double wv2 = waveVectors[i].squared();
                for (int j = 0; j < wv2List.size(); j++) {
                    if (Math.abs(wv2List.get(j) - wv2) < 1e-6) {
                        map[i] = j;
                        break;
                    }
                }
                if (map[i] >= 0) continue;
                map[i] = wv2List.size();
                wv2List.add(wv2);
            } else {
                double wv1 = Math.abs(waveVectors[i].getX(motionDim));
                double wv22 = 0;
                for (int j = 0; j < waveVectors[i].getD(); j++) {
                    if (j == motionDim) continue;
                    double wvj = waveVectors[i].getX(j);
                    wv22 += wvj * wvj;
                }
                for (int j = 0; j < wv2List.size(); j++) {
                    if (Math.abs(wv2List.get(j) - wv22) < 1e-6 && Math.abs(wv1List.get(j) - wv1) < 1e-6) {
                        map[i] = j;
                        break;
                    }
                }
                if (map[i] >= 0) continue;
                map[i] = wv2List.size();
                wv2List.add(wv22);
                wv1List.add(wv1);
            }
        }
        return map;
    }

    public void setMinInterval(int minInterval) {
        this.minInterval = minInterval;
    }

    public void reset() {
        reallocate(0);
    }

    protected void reallocate(int n) {
        int nOld = skipNow.length;
        skipNow = Arrays.copyOf(skipNow, n);
        tData = new DataDoubleArray(new int[]{n});
        tDataInfo = new DataDoubleArray.DataInfoDoubleArray("t", Time.DIMENSION, new int[]{n});
        tDataInfo.addTag(tTag);
        for (int i = 0; i < corSum.length; i++) {
            corSum[i] = Arrays.copyOf(corSum[i], n);
            sum2[i] = Arrays.copyOf(sum2[i], n);
            nSamplesCor[i] = Arrays.copyOf(nSamplesCor[i], n);
            nSamples2[i] = Arrays.copyOf(nSamples2[i], n);
            data[i] = new DataFunction(new int[]{n});
            dataInfo[i] = new DataFunction.DataInfoFunction("alpha", Null.DIMENSION, this);
            dataInfo[i].addTag(tag[i]);
        }
        for (int i = 0; i < lastXY.length; i++) {
            lastXY[i] = Arrays.copyOf(lastXY[i], n);
            for (int j = nOld; j < n; j++) {
                lastXY[i][j] = new double[]{Double.NaN, Double.NaN};
                skipNow[j] = true;
            }
        }
        setTimeData();
    }

    protected void setTimeData() {
        double[] t = tData.getData();
        if (t.length > 0) {
            double[] savedTimes = configStorage.getSavedTimes();
            if (savedTimes[1] < 0) {
                Arrays.fill(t, Double.NaN);
            } else {
                double dt = savedTimes[0] - savedTimes[1];
                for (int i = 0; i < t.length; i++) {
                    t[i] = dt * (1L << i);
                }
            }
        }
    }

    public void putData(int j, double[][] xyData) {
        if (Double.isNaN(xyData[0][0])) return;
        if (j >= corSum[0].length) reallocate(j + 1);
        else if (Double.isNaN(tData.getValue(0))) {
            setTimeData();
        }
        // we receive (x^2+y^2)
        for (int i = 0; i < xyData.length; i++) {
            double x = xyData[i][0];
            double y = xyData[i][1];
            int imap = wvMap[i];

            if (!skipNow[j]) {
                corSum[imap][j] += x * lastXY[i][j][0] + y * lastXY[i][j][1];
                nSamplesCor[imap][j]++;
            }
            sum2[imap][j] += x * x + y * y;
            nSamples2[imap][j]++;
            lastXY[i][j][0] = x;
            lastXY[i][j][1] = y;
        }
        if (skipNow[j]) {
            skipNow[j] = false;
        } else if (j < minInterval) {
            skipNow[j] = true;
        }
    }

    public void putData(int j, IData inputData, double[] phaseAngles) {
        if (Double.isNaN(inputData.getValue(0))) return;
        double[][] xyData = new double[phaseAngles.length][];
        for (int i = 0; i < data.length; i++) {
            double sfac = inputData.getValue(i);
            double tanphi = Math.tan(phaseAngles[i]);
            double x = Math.sqrt(sfac / (1 + tanphi * tanphi));
            if (phaseAngles[i] > Math.PI / 2 || phaseAngles[i] < -Math.PI / 2) {
                x = -x;
            }
            double y = x * tanphi;
            xyData[i] = new double[]{x, y};
        }
        putData(j, xyData);
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

    public Sink makeSink(int idx, MeterStructureFactor meterSFac) {
        return new Sink(idx, meterSFac);
    }

    public Meter makeMeter(int idx) {
        return new Meter(idx);
    }

    public class Sink implements IDataSink {

        private final int idx;
        private final MeterStructureFactor meterSFac;

        public Sink(int idx, MeterStructureFactor meterSFac) {
            this.idx = idx;
            this.meterSFac = meterSFac;
        }

        @Override
        public void putData(IData data) {
            StructureFactorComponentCorrelation.this.putData(idx, data, meterSFac.getPhaseAngles());
        }

        public void putDataInfo(IDataInfo inputDataInfo) {
            if (inputDataInfo.getLength() != data.length) {
                throw new RuntimeException("# of wave vectors changed!");
            }
        }
    }

    public class Meter implements IDataSource {

        private final int idx;

        public Meter(int idx) {
            this.idx = idx;
        }

        @Override
        public IData getData() {
            double[] y = data[idx].getData();
            for (int i = 0; i < y.length; i++) {
                if (nSamplesCor[idx][i] < 1) {
                    y[i] = Double.NaN;
                    continue;
                }
                y[i] = corSum[idx][i] / nSamplesCor[idx][i] / (sum2[idx][i] / nSamples2[idx][i]);
            }
            return data[idx];
        }

        @Override
        public DataTag getTag() {
            return tag[idx];
        }

        @Override
        public IDataInfo getDataInfo() {
            return dataInfo[idx];
        }
    }
}
