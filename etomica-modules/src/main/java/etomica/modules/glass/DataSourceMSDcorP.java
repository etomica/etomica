package etomica.modules.glass;

import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.integrator.IntegratorMD;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Time;

import java.util.Arrays;

public class DataSourceMSDcorP implements IDataSink, IDataSource, DataSourceIndependent, DataSourceMSD.MSDSink {

    protected DataDoubleArray tData;
    protected DataDoubleArray.DataInfoDoubleArray tDataInfo;
    protected DataFunction data;
    protected DataFunction.DataInfoFunction dataInfo;
    protected double[] msdSum, pSum, pmsdSum, msd2Sum, p2Sum;
    protected final DataTag tTag, tag;
    protected long[] nSamples;
    protected double[] lastSampleMSD, lastSampleP;
    protected double[] blockSumP;
    protected long[] nBlockSamplesP;
    protected long[] lastStepMSD, lastStepP;
    protected final IntegratorMD integrator;
    protected long step0;
    protected boolean enabled;

    public DataSourceMSDcorP(IntegratorMD integrator) {
        this.integrator = integrator;
        lastSampleMSD = lastSampleP = pSum = pmsdSum = msdSum = msd2Sum = p2Sum = new double[0];
        nBlockSamplesP = new long[60];
        lastStepP = new long[60];
        blockSumP = new double[60];
        lastStepMSD = nSamples = new long[0];
        tag = new DataTag();
        tTag = new DataTag();
        resetStep0();
    }

    public void resetStep0() {
        step0 = integrator.getStepCount();
        for (int i = 0; i < nBlockSamplesP.length; i++) {
            lastStepP[i] = nBlockSamplesP[i] = 0;
            blockSumP[i] = 0;
        }
        reallocate(0);
    }

    public void reallocate(int n) {
        if (n == msdSum.length && data != null) return;
        msdSum = Arrays.copyOf(msdSum, n);
        msd2Sum = Arrays.copyOf(msd2Sum, n);
        pSum = Arrays.copyOf(pSum, n);
        p2Sum = Arrays.copyOf(p2Sum, n);
        pmsdSum = Arrays.copyOf(pmsdSum, n);
        nSamples = Arrays.copyOf(nSamples, n);
        lastStepMSD = Arrays.copyOf(lastStepMSD, n);
        lastSampleP = Arrays.copyOf(lastSampleP, n);
        lastSampleMSD = Arrays.copyOf(lastSampleMSD, n);
        data = new DataFunction(new int[]{n});
        tData = new DataDoubleArray(new int[]{n});
        tDataInfo = new DataDoubleArray.DataInfoDoubleArray("t", Time.DIMENSION, new int[]{n});
        tDataInfo.addTag(tTag);
        dataInfo = new DataFunction.DataInfoFunction("P-MSD cor", Null.DIMENSION, this);
        dataInfo.addTag(tag);
        double[] t = tData.getData();
        if (t.length > 0) {
            double dt = integrator.getTimeStep();
            for (int i = 0; i < t.length; i++) {
                t[i] = dt * (1L << i);
            }
        }
    }

    public void setEnabled(boolean isEnabled) {
        enabled = isEnabled;
        resetStep0();
    }

    public boolean getEnabled() {
        return enabled;
    }

    @Override
    public IData getData() {
        double[] y = data.getData();
        for (int i = 0; i < y.length - 1; i++) {
            double pAvg = pSum[i] / nSamples[i];
            double msdAvg = msdSum[i] / nSamples[i];
            double pmsdAvg = pmsdSum[i] / nSamples[i];
            double p2Avg = p2Sum[i] / nSamples[i];
            double pVar = p2Avg - pAvg * pAvg;
            double msd2Avg = msd2Sum[i] / nSamples[i];
            double msdVar = msd2Avg - msdAvg * msdAvg;
            y[i] = (pmsdAvg - pAvg * msdAvg) / Math.sqrt(msdVar * pVar);
        }
        if (y.length > 0) y[y.length - 1] = Double.NaN;
        return data;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    // pressure
    @Override
    public void putData(IData inputData) {
        if (!enabled) return;
        double p = inputData.getValue(0);
        long step = integrator.getStepCount() - step0;
        for (int i = 0; i < blockSumP.length; i++) {
            blockSumP[i] += p;
            nBlockSamplesP[i]++;
            if (step % (1L << i) == 0) {
                if (lastSampleP.length <= i) reallocate(i + 1);
                lastSampleP[i] = blockSumP[i] / nBlockSamplesP[i];
                lastStepP[i] = step;
                blockSumP[i] = 0;
                nBlockSamplesP[i] = 0;
                if (lastStepMSD[i] == step) {
                    processSample(i);
                }
            }
        }
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

    public void putMSD(int log2interval, long step, double msd) {
        if (lastSampleMSD.length <= log2interval) reallocate(log2interval + 1);
        lastSampleMSD[log2interval] = msd;
        if (lastStepP[log2interval] == step) {
            processSample(log2interval);
        }
    }

    private void processSample(int log2interval) {
        pSum[log2interval] += lastSampleP[log2interval];
        msdSum[log2interval] += lastSampleMSD[log2interval];
        pmsdSum[log2interval] += lastSampleP[log2interval] * lastSampleMSD[log2interval];
        p2Sum[log2interval] += lastSampleP[log2interval] * lastSampleP[log2interval];
        msd2Sum[log2interval] += lastSampleMSD[log2interval] * lastSampleMSD[log2interval];
        lastStepMSD[log2interval] = lastStepP[log2interval] = -1;
        nSamples[log2interval]++;
    }
}
