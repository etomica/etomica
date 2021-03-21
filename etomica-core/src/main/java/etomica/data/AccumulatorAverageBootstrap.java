/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.histogram.Histogram;
import etomica.data.histogram.HistogramCollapsing;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.math.function.Function;
import etomica.math.function.FunctionDifferentiable;
import etomica.math.function.FunctionInvertible;
import etomica.math.function.IFunction;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Quantity;
import etomica.util.random.IRandom;
import etomica.util.random.RandomMersenneTwister;
import etomica.util.random.RandomNumberGenerator;
import etomica.util.random.RandomNumberGeneratorUnix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


/**
 * AccumulatorAverage that keeps track of block averages and standard
 * deviations corresponding to all block sizes from 1 to the n where n is the
 * largest power of 2 less than the current number of samples.  In addition to
 * block averages and standard deviations of the sampled quantity, a function
 * can be given to the class such that it will compute pass block averages
 * through that function and compute standard deviation of the transformed
 * block averages (this is what makes this class actually useful).  The
 * transformed values can be obtained via conventiently named methods with
 * "Log" appended to end.  This class can also compute histograms of the block
 * averages for all the different block sizes.  The histograms are not returned
 * by any method, but are instead pushed to data sinks.
 * <p>
 * This class also uses a method to improve the block statistics for large
 * block sizes based on bootstrapping.  All block averages are kept for a
 * particular block size (n1 blocks of size b1).  Statistics for that block
 * are computed from that data.  Block averages for blocks of the next larger
 * size (2*n1) are computed by randomly selecting n1/2 pairs of block values
 * from the b1 set (without including any block value more than once).
 * After using all block values, n1/2 additional pairs are taken from the
 * full b1 set, such that n1 blocks of size b1*2 are obtained, instead of
 * n1/2 blocks that would normally exist.  Statistics for block size b1*2 are
 * computed from these n1 blocks.  This process is repeated for each block size
 * such that the largest block size (which would normally have 1 block) has n1
 * blocks.
 * <p>
 * The class can operate with or without replacement during replacement.
 * Replacement is a standard part of bootstrapping, but does not seem to have a
 * significant impact here due to the way data is built up from small to large
 * block sizes.  If replacement is turned off, the resampling goes
 * substantially faster.
 * <p>
 * This accumulator can only operate on Data with a single value.
 */
public class AccumulatorAverageBootstrap extends DataAccumulator implements DataSourceIndependent {

    protected final DataTag nTag;
    protected final IRandom random;
    protected final int nMoments;
    protected boolean withReplacement = true;
    protected long count;
    protected double[] blockSums;
    protected double[][] sums, lSums;
    protected DataFunction averages, stdev, skew;
    protected DataInfoDoubleArray nDataInfo;
    protected DataDoubleArray nData;
    protected double[] rawData, rawData2, rawData3;
    protected int nRawData, rawDataBlockSize;
    protected int nRawDataDoubles, maxNumBlocks;
    protected List<Integer> intArrayList;
    protected int initialSeed;
    protected AccumulatorAverageCollapsing lAvg;
    protected DataDouble lData;
    protected DataInfoDouble lDataInfo;
    protected List<Histogram> histogramList;
    protected HistogramFactory histogramFactory;
    protected final List<DataFunction> histograms;
    protected final List<DataInfoFunction> histogramsInfo;
    protected final ArrayList<IDataSink> histogramSinks;
    protected IFunction transformFunction;
    protected boolean doUntransform;

    public AccumulatorAverageBootstrap() {
        this(2);
    }

    /**
     * With this constructor, the accumulator will use an internal RNG with a
     * set seed such that each resampling of the data will be identical (the
     * same set of samples will go into the averages and the standard
     * deviations, etc).
     *
     * @param nMoments the number of moments handled by the accumulator.
     *                 should be 2 (standard deviation) or 3 (stdev + skew).
     */
    public AccumulatorAverageBootstrap(int nMoments) {
        this(new RandomMersenneTwister(RandomNumberGeneratorUnix.getRandSeedArray()), nMoments);
        initialSeed = random.nextInt(1 << 30);
    }

    public AccumulatorAverageBootstrap(IRandom random) {
        this(random, 2);
    }

    /**
     * @param nMoments the number of moments handled by the accumulator.
     *                 should be 2 (standard deviation) or 3 (stdev + skew).
     */
    public AccumulatorAverageBootstrap(IRandom random, int nMoments) {
        if (nMoments < 2 || nMoments > 3) {
            throw new IllegalArgumentException("nMoments should be 2 or 3");
        }
        this.nMoments = nMoments;
        this.random = random;
        initialSeed = -1;
        setPushInterval(100);
        setNumRawDataDoubles(20);
        setMaxNumBlocks(12);
        nTag = new DataTag();
        histogramList = new ArrayList<>();
        histogramFactory = new HistogramFactory() {
            @Override
            public Histogram makeHistogram() {
                return new HistogramCollapsing();
            }
        };
        transformFunction = Function.Identity.INSTANCE;
        histograms = new ArrayList<>();
        histogramsInfo = new ArrayList<>();
        histogramSinks = new ArrayList<>();
    }

    /**
     * @param newFunction the transform function, which block averages are passed through.
     */
    public void setTransformFunction(IFunction newFunction) {
        transformFunction = newFunction;
    }

    /**
     * @param doUntransform causes histograms to be untransformed once they are constructed
     */
    public void setUntransformHistograms(boolean doUntransform) {
        this.doUntransform = doUntransform;
    }

    public static interface HistogramFactory {
        public Histogram makeHistogram();
    }

    /**
     * @param newFactory the factory that makes histograms
     */
    public void setHistogramFactory(HistogramFactory newFactory) {
        this.histogramFactory = newFactory;
    }

    /**
     * Sets the sink that will receive histograms for a block when pushHistograms is called.
     */
    public void setHistogramSink(int blockSize, IDataSink sink) {
        while (histogramSinks.size() <= blockSize) histogramSinks.add(null);
        histogramSinks.set(blockSize, sink);
    }

    /**
     * @param newWithReplacement enables sampling with replacement
     */
    public void setWithReplacement(boolean newWithReplacement) {
        withReplacement = newWithReplacement;
    }

    public static void main(String[] args) {
        AccumulatorAverageBootstrap ac = new AccumulatorAverageBootstrap(new RandomNumberGenerator());
        ac.putDataInfo(new DataInfoDouble("foo", Null.DIMENSION));
        DataDouble d = new DataDouble();
        RandomNumberGenerator rng = new RandomNumberGenerator();
        for (long i = 0; i < 1L << 25; i++) {
            d.x = Math.exp(rng.nextGaussian());
            ac.putData(d);
        }
        System.out.println(ac.getAverages().getValue(0) + " " + ac.getAverages().getValue(1) + " " + ac.getAverages().getValue(2) + " " + ac.getAverages().getValue(3));
        System.out.println(ac.getStdev().getValue(0) + " " + ac.getStdev().getValue(1) + " " + ac.getStdev().getValue(2) + " " + ac.getStdev().getValue(3));
        System.out.println(ac.getAverageLogs().getValue(0) + " " + ac.getAverageLogs().getValue(1) + " " + ac.getAverageLogs().getValue(2) + " " + ac.getAverageLogs().getValue(3));
        System.out.println(ac.getStdevLog().getValue(0) + " " + ac.getStdevLog().getValue(1) + " " + ac.getStdevLog().getValue(2) + " " + ac.getStdevLog().getValue(3));
    }

    /**
     * Sets the number of times the random block generation method is used.
     * The number of blocks used will be 2^numRawDataDoubles.
     */
    public void setNumRawDataDoubles(int newNumRawDataDoubles) {
        nRawDataDoubles = newNumRawDataDoubles;
        intArrayList = new ArrayList<Integer>(1 << nRawDataDoubles);
        reset();
    }

    /**
     * Sets the maximum number of blocks of raw data to be 1<<maxNumBlocks
     */
    public void setMaxNumBlocks(int newMaxNumBlocks) {
        maxNumBlocks = newMaxNumBlocks;
        reset();
    }

    public IData processData(IData inputData) {
        if (!active) return null;
        int oldSize = sums.length;
        addData(inputData);
        if (sums.length > oldSize && false) {
            return getData();
        }
        return null;
    }

    /**
     * Add the given values to the sums and block sums. If any of the given data
     * values is NaN, method returns with no effect on accumulation sums.
     */
    public boolean addData(IData data) {
        if (data.isNaN())
            return false;
        double value = data.getValue(0);
        if (value == Double.POSITIVE_INFINITY) {
            throw new RuntimeException("got " + value);
        }
        if (value == 0) {
            value = Double.MIN_VALUE;
        }
        if (lAvg != null) {
            lData.x = transformFunction.f(value);
            lAvg.putData(lData);
        }
        if (rawDataBlockSize == 0) {
            // special case early part where we keep every sample as raw data
            rawData[nRawData] = value;
            if (histogramList.size() == 0) {
                histogramList.add(histogramFactory.makeHistogram());
            }
            histogramList.get(0).addValue(transformFunction.f(value));
            nRawData++;
            count++;
            int nDoubles = 63 - Long.numberOfLeadingZeros(count);
            if (nDoubles >= averages.getLength()) {
                resizeData(nDoubles + 1);
            }
            if (nRawData == 1 << maxNumBlocks) {
                // we filled our raw data.  do summation for blockSize=0,
                // prepare for blockSize=1
                histogramList.add(histogramFactory.makeHistogram());
                for (int i = 0; i < nRawData; i += 2) {
                    doSums(sums[0], rawData[i]);
                    doSums(lSums[0], transformFunction.f(rawData[i]));
                    doSums(sums[0], rawData[i + 1]);
                    doSums(lSums[0], transformFunction.f(rawData[i + 1]));
                    rawData[i / 2] = 0.5 * (rawData[i] + rawData[i + 1]);
                    histogramList.get(1).addValue(transformFunction.f(rawData[i / 2]));
                }
                for (int i = nRawData / 2; i < nRawData; i++) {
                    rawData[i] = 0;
                }
                rawDataBlockSize++;
                nRawData /= 2;
            }
            return true;
        }
        // process sample as blockSize=0 data
        doSums(sums[0], value);
        doSums(lSums[0], transformFunction.f(value));
        histogramList.get(0).addValue(transformFunction.f(value));
        blockSums[1] += value;
        count++;
        int nDoubles = Long.numberOfTrailingZeros(count);
        for (int i = 1; i < nDoubles + 1; i++) {
            blockSums[i] *= 0.5;
            if (rawDataBlockSize == i) {
                // we consider this block average to be raw data
                rawData[nRawData] = blockSums[i];
                blockSums[i] = 0;
                nRawData++;
                if (nRawData == 1 << maxNumBlocks) {
                    // we filled our raw data
                    histogramList.add(histogramFactory.makeHistogram());
                    resizeData(nDoubles + 1);
                    for (int j = 0; j < nRawData; j += 2) {
                        doSums(sums[i], rawData[j]);
                        doSums(lSums[i], transformFunction.f(rawData[j]));
                        doSums(sums[i], rawData[j + 1]);
                        doSums(lSums[i], transformFunction.f(rawData[j + 1]));
                        rawData[j / 2] = 0.5 * (rawData[j] + rawData[j + 1]);
                        histogramList.get(i + 1).addValue(transformFunction.f(rawData[j / 2]));
                    }
                    for (int j = nRawData / 2; j < nRawData; j++) {
                        rawData[j] = 0;
                    }
                    rawDataBlockSize++;
                    nRawData /= 2;
                }
                return true;
            }
            // process block average normally
            doSums(sums[i], blockSums[i]);
            doSums(lSums[i], transformFunction.f(blockSums[i]));
            blockSums[i + 1] += blockSums[i];
            blockSums[i] = 0;
        }
        return true;
    }

    public long getCount() {
        return count;
    }

    protected void resizeData(int newSize) {
        blockSums = Arrays.copyOf(blockSums, newSize + 1);
        sums = Arrays.copyOf(sums, newSize);
        sums[newSize - 1] = new double[nMoments];
        lSums = Arrays.copyOf(lSums, newSize);
        lSums[newSize - 1] = new double[nMoments];
        averages = new DataFunction(new int[]{newSize});
        stdev = new DataFunction(new int[]{newSize});
        if (nMoments > 2) {
            skew = new DataFunction(new int[]{newSize});
        }
        for (int i = histograms.size(); i < newSize; i++) {
            histogramList.add(histogramFactory.makeHistogram());
        }
        dataInfo = new DataInfoFunction(dataInfo.getLabel(), dataInfo.getDimension(), this);
        nDataInfo = new DataInfoDoubleArray("block size", Quantity.DIMENSION, new int[]{newSize});
        nDataInfo.addTag(nTag);
        nData = new DataDoubleArray(newSize);
        double[] n = nData.getData();
        for (int i = 0; i < n.length; i++) {
            n[i] = 1L << i;
        }

        if (dataSink != null) {
            dataSink.putDataInfo(dataInfo);
        }
    }

    /*
     * Returns the average logs
     */
    public IData getData() {
        return getAverageLogs();
    }

    protected void untransform(double[] x, double[] h, IFunction f) {
        int n = x == null ? h.length : x.length;
        for (int i = 0; i < n; i++) {
            if (x != null) x[i] = ((FunctionInvertible) f).inverse(x[i]);
            if (h != null) h[i] = h[i] * ((FunctionDifferentiable) f).df(1, x[i]);
        }
    }

    protected void initHistogram(int blockSize, IDataSink sink) {
        DataFunction histogramData = new DataFunction(new int[]{0}, new double[0]);
        histograms.set(blockSize, histogramData);
        double[] oldX = new double[0];
        DataInfoDoubleArray xDataInfo = new DataInfoDoubleArray("x bs " + (1L << blockSize), Null.DIMENSION, new int[]{0});
        DataSourceIndependentSimple xDataSource = new DataSourceIndependentSimple(oldX, xDataInfo);
        DataInfoFunction histogramDataInfo = new DataInfoFunction("bs " + (1L << blockSize), Null.DIMENSION, xDataSource);
        // XXX we need to add incoming tags
        // XXX we should have a separate tag (additional) for each block size
        histogramDataInfo.addTag(tag);
        histogramsInfo.set(blockSize, histogramDataInfo);
        if (sink != null) {
            sink.putDataInfo(histogramDataInfo);
        }
    }

    /**
     * Constructs and pushes histograms to any histogram data sink.
     */
    public void pushHistograms() {
        for (int blockSize = 0; blockSize <= rawDataBlockSize; blockSize++) {
            if (histogramSinks.size() < blockSize + 1) return;
            IDataSink sink = histogramSinks.get(blockSize);
            while (histograms.size() <= blockSize) {
                histograms.add(null);
                histogramsInfo.add(null);
            }
            DataFunction histogramData = histograms.get(blockSize);
            DataInfoFunction histogramDataInfo = histogramsInfo.get(blockSize);
            if (histogramData == null) {
                initHistogram(blockSize, sink);
                histogramData = histograms.get(blockSize);
                histogramDataInfo = histogramsInfo.get(blockSize);
            }
            double[] x = histogramList.get(blockSize).xValues();
            double[] oldX = histogramDataInfo.getXDataSource().getIndependentData(0).getData();
            if (oldX.length != x.length) {
                oldX = x.clone();
                DataInfoDoubleArray xDataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{oldX.length});
                ((DataSourceIndependentSimple) histogramDataInfo.getXDataSource()).update(oldX, xDataInfo);
            } else {
                System.arraycopy(x, 0, oldX, 0, x.length);
            }
            double[] h = histogramList.get(blockSize).getHistogram();
            double[] oldH = histogramData.getData();
            if (h.length != oldH.length) {
                oldH = h.clone();
                histogramData = new DataFunction(new int[]{oldH.length}, oldH);
            } else {
                System.arraycopy(h, 0, oldH, 0, h.length);
            }
            if (doUntransform) {
                untransform(oldX, oldH, transformFunction);
            }
            if (sink != null) {
                sink.putData(histogramData);
            }
        }

        for (int j = 0; j < nRawData; j++) {
            rawData2[j] = rawData[j];
        }

        for (int i = rawDataBlockSize + 1; i < rawDataBlockSize + nRawDataDoubles; i++) {
            if (histogramSinks.size() < i + 1) return;
            reblockData();
            while (histograms.size() <= i) {
                histograms.add(null);
                histogramsInfo.add(null);
            }
            Histogram hist = histogramList.get(i);
            hist.reset();
            for (int j = 0; j < nRawData; j++) {
                hist.addValue(transformFunction.f(rawData2[j]));
            }

            IDataSink sink = histogramSinks.get(i);
            DataFunction histogramData = histograms.get(i);
            DataInfoFunction histogramDataInfo = histogramsInfo.get(i);
            if (histogramData == null) {
                initHistogram(i, sink);
                histogramData = histograms.get(i);
                histogramDataInfo = histogramsInfo.get(i);
            }
            double[] oldX = histogramDataInfo.getXDataSource().getIndependentData(0).getData();
            double[] x = histogramList.get(i).xValues();
            if (oldX.length != x.length) {
                oldX = x.clone();
                DataInfoDoubleArray xDataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{oldX.length});
                ((DataSourceIndependentSimple) histogramDataInfo.getXDataSource()).update(oldX, xDataInfo);
            } else {
                System.arraycopy(x, 0, oldX, 0, x.length);
            }

            double[] h = histogramList.get(i).getHistogram();
            double[] oldH = histogramData.getData();
            if (h.length != oldH.length) {
                oldH = h.clone();
                histogramData = new DataFunction(new int[]{oldH.length}, oldH);
            }
            if (doUntransform) {
                untransform(oldX, oldH, transformFunction);
            }
            if (sink != null) {
                sink.putData(histogramData);
            }
        }
    }

    /**
     * Returns the averages for all block sizes
     */
    public IData getAverages() {
        double[] av = averages.getData();
        for (int i = 0; i < rawDataBlockSize; i++) {
            av[i] = sums[i][0] / (count >> i);
        }
        int n = averages.getLength();
        double iSum = 0;
        for (int j = 0; j < nRawData; j++) {
            iSum += rawData[j];
        }
        for (int i = rawDataBlockSize; i < n; i++) {
            av[i] = iSum / nRawData;
        }
        return averages;
    }

    /**
     * Returns the average natural logs of the block averages, for all block
     * sizes.
     */
    public IData getAverageLogs() {
        double[] av = averages.getData();
        if (av.length == 0) {
            return averages;
        }
        if (initialSeed >= 0) {
            ((RandomMersenneTwister) random).setSeed(initialSeed);
        }
        for (int i = 0; i < rawDataBlockSize; i++) {
            av[i] = lSums[i][0] / (count >> i);
        }
        double lSum = 0;
        for (int j = 0; j < nRawData; j++) {
            lSum += transformFunction.f(rawData[j]);
            rawData2[j] = rawData[j];
        }
        av[rawDataBlockSize] = lSum / nRawData;
        for (int i = rawDataBlockSize + 1; i < rawDataBlockSize + (32 - Integer.numberOfLeadingZeros(nRawData)); i++) {
            lSum = 0;
            reblockData();
            for (int j = 0; j < nRawData; j++) {
                lSum += transformFunction.f(rawData2[j]);
            }
            av[i] = lSum / nRawData;
        }
        return averages;
    }

    /**
     * Reblock data from rawData2 into rawData3, then copy it back to rawData2
     */
    protected void reblockData() {
        if (!withReplacement) {
            // remove elements as we select them
            intArrayList.clear();
            for (int j = 0; j < nRawData; j++) {
                ensureHappyIntArrayList();
                int idxidx1 = random.nextInt(intArrayList.size());
                int idx1 = intArrayList.get(idxidx1);
                intArrayList.remove(idxidx1);
                ensureHappyIntArrayList();
                int idxidx2 = random.nextInt(intArrayList.size());
                int idx2 = intArrayList.get(idxidx2);
                intArrayList.remove(idxidx2);
                double sum12 = (rawData2[idx1] + rawData2[idx2]) / 2;
                rawData3[j] = sum12;
            }
        } else {
            // with replacement
            ensureHappyIntArrayList();
            for (int j = 0; j < nRawData; j++) {
                int idx1 = random.nextInt(nRawData);
                int idx2 = random.nextInt(nRawData);
                rawData3[j] = 0.5 * (rawData2[idx1] + rawData2[idx2]);
            }
        }
        for (int j = 0; j < nRawData; j++) {
            rawData2[j] = rawData3[j];
        }
    }

    /**
     * Brute-force reblock data from rawData into rawData2,
     * taking blocks of size n
     */
    protected void reblockData(int n) {
        // with replacement
        for (int j = 0; j < nRawData; j++) {
            double innerSum = 0;
            for (int i = 0; i < n; i++) {
                innerSum += rawData[random.nextInt(nRawData)];
            }
            innerSum /= n;
            rawData2[j] = innerSum;
        }
    }

    /**
     * Refills inArrayList if it's empty.
     */
    protected void ensureHappyIntArrayList() {
        if (intArrayList.size() == 0) {
            for (int j = 0; j < nRawData; j++) {
                intArrayList.add(j);
            }
        }
    }

    /**
     * Returns the standard deviation of the block averages, for all block
     * sizes.
     */
    public IData getStdev() {
        if (initialSeed >= 0) {
            ((RandomMersenneTwister) random).setSeed(initialSeed);
        }
        double[] std = stdev.getData();
        if (std.length == 0) {
            return stdev;
        }
        for (int i = 0; i < rawDataBlockSize; i++) {
            long iCount = count >> i;
            double av = sums[i][0] / iCount;
            std[i] = Math.sqrt(sums[i][1] / iCount - av * av);
            if (Double.isNaN(std[i])) {
                std[i] = 0;
            }
        }
        double iSum = 0, iSum2 = 0;
        for (int j = 0; j < nRawData; j++) {
            iSum += rawData[j];
            iSum2 += rawData[j] * rawData[j];
            rawData2[j] = rawData[j];
        }
        iSum /= nRawData;
        iSum2 /= nRawData;
        double diff = iSum2 - iSum * iSum;
        if (diff < 0) {
            std[rawDataBlockSize] = 0;
        } else {
            std[rawDataBlockSize] = Math.sqrt(diff);
        }
        for (int i = rawDataBlockSize + 1; i < rawDataBlockSize + (32 - Integer.numberOfLeadingZeros(nRawData)); i++) {
            iSum = 0;
            iSum2 = 0;
            reblockData();
            for (int j = 0; j < nRawData; j++) {
                iSum += rawData2[j];
                iSum2 += rawData2[j] * rawData2[j];
            }
            iSum /= nRawData;
            iSum2 /= nRawData;
            diff = iSum2 - iSum * iSum;
            if (diff < 0) {
                std[rawDataBlockSize] = 0;
            } else {
                std[i] = Math.sqrt(diff);
            }
        }
        return stdev;
    }

    /**
     * Returns the standard deviation of the log of the block averages, for all
     * block sizes.
     */
    public IData getStdevLog() {
        if (initialSeed >= 0) {
            ((RandomMersenneTwister) random).setSeed(initialSeed);
        }
        double[] std = stdev.getData();
        if (std.length == 0) {
            return stdev;
        }
        for (int i = 0; i < rawDataBlockSize; i++) {
            long iCount = count >> i;
            double av = lSums[i][0] / iCount;
            std[i] = Math.sqrt(lSums[i][1] / iCount - av * av);
        }
        double iSum = 0, iSum2 = 0;
        for (int j = 0; j < nRawData; j++) {
            double l = transformFunction.f(rawData[j]);
            iSum += l;
            iSum2 += l * l;
            rawData2[j] = rawData[j];
        }
        std[rawDataBlockSize] = Math.sqrt(iSum2 / nRawData - (iSum * iSum) / (nRawData * nRawData));
        for (int i = rawDataBlockSize + 1; i < rawDataBlockSize + (32 - Integer.numberOfLeadingZeros(nRawData)); i++) {
            iSum = 0;
            iSum2 = 0;
            reblockData();
            for (int j = 0; j < nRawData; j++) {
                double l = transformFunction.f(rawData2[j]);
                iSum += l;
                iSum2 += l * l;
            }
            std[i] = Math.sqrt(iSum2 / nRawData - (iSum * iSum) / (nRawData * nRawData));
        }
        return stdev;
    }

    /**
     * Returns the standard deviation of the log of the block averages, for all
     * block sizes.
     */
    public IData getSkewLog() {
        if (initialSeed >= 0) {
            ((RandomMersenneTwister) random).setSeed(initialSeed);
        }
        if (nMoments < 3) {
            throw new RuntimeException("skew not handled");
        }
        double[] std = skew.getData();
        for (int i = 0; i < rawDataBlockSize; i++) {
            long iCount = count >> i;
            double av = lSums[i][0] / iCount;
            double av2 = lSums[i][1] / iCount;
            double std3 = Math.pow(av2 - av * av, 1.5);
            std[i] = (lSums[i][2] / iCount - 3 * av * av2 + 2 * av * av * av) / std3;
        }
        double iSum = 0, iSum2 = 0, iSum3 = 0;
        for (int j = 0; j < nRawData; j++) {
            double l = transformFunction.f(rawData[j]);
            iSum += l;
            iSum2 += l * l;
            iSum3 += l * l * l;
            rawData2[j] = rawData[j];
        }
        double av = iSum / nRawData;
        double av2 = iSum2 / nRawData;
        double std3 = Math.pow(av2 - av * av, 1.5);
        std[rawDataBlockSize] = (iSum3 / nRawData - 3 * av * av2 + 2 * av * av * av) / std3;
        int iMax = rawDataBlockSize + (32 - Integer.numberOfLeadingZeros(nRawData));
        for (int i = rawDataBlockSize + 1; i < iMax; i++) {
            iSum = 0;
            iSum2 = 0;
            iSum3 = 0;
            reblockData();
            for (int j = 0; j < nRawData; j++) {
                double l = transformFunction.f(rawData2[j]);
                iSum += l;
                iSum2 += l * l;
                iSum3 += l * l * l;
            }
            av = iSum / nRawData;
            av2 = iSum2 / nRawData;
            std3 = Math.pow(av2 - av * av, 1.5);
            std[i] = (iSum3 / nRawData - 3 * av * av2 + 2 * av * av * av) / std3;
        }
        return skew;
    }

    /**
     * Performs the block sum after
     */
    protected void doSums(double[] subSums, double v) {
        double vpj = 1;
        for (int i = 0; i < nMoments; i++) {
            vpj *= v;
            subSums[i] += vpj;
        }
    }

    public int getNumRawData() {
        return nRawData;
    }

    public double[] getRawData() {
        return rawData;
    }

    public double getSampleVariance() {
        if (rawDataBlockSize > 0) {
            return (count * sums[0][1] - sums[0][0] * sums[0][0]) / (count * (count - 1));
        }
        double sums0 = 0, sums1 = 0;
        for (int i = 0; i < nRawData; i++) {
            sums0 += rawData[i];
            sums1 += rawData[i] * rawData[i];
        }
        return (count * sums1 - sums0 * sums0) / (count * (count - 1));
    }

    public double getSampleLogVariance() {
        if (rawDataBlockSize > 0) {
            return (count * lSums[0][1] - lSums[0][0] * lSums[0][0]) / (count * (count - 1));
        }
        double sums0 = 0, sums1 = 0;
        for (int i = 0; i < nRawData; i++) {
            double l = transformFunction.f(rawData[i]);
            sums0 += l;
            sums1 += l * l;
        }
        return (count * sums1 - sums0 * sums0) / (count * (count - 1));
    }

    public int getRawDataBlockSize() {
        return rawDataBlockSize;
    }

    /**
     * Resets all sums to zero. All statistics are cleared.
     */
    public void reset() {
        count = 0;
        blockSums = new double[1];
        sums = new double[0][nMoments];
        lSums = new double[0][nMoments];
        averages = new DataFunction(new int[]{0});
        stdev = new DataFunction(new int[]{0});
        if (nMoments > 2) {
            skew = new DataFunction(new int[]{0});
        }
        if (dataInfo != null) {
            dataInfo = new DataInfoFunction(dataInfo.getLabel(), dataInfo.getDimension(), this);
        }
        nDataInfo = new DataInfoDoubleArray("block size", Quantity.DIMENSION, new int[]{0});
        nDataInfo.addTag(nTag);
        rawData = new double[1 << maxNumBlocks];
        rawData2 = new double[1 << maxNumBlocks];
        rawData3 = new double[1 << maxNumBlocks];
        rawDataBlockSize = 0;
        nRawData = 0;
        if (histogramSinks != null) {
            for (int i = 0; i < histogramSinks.size(); i++) {
                initHistogram(i, histogramSinks.get(i));
                if (histogramSinks.get(i) != null) {
                    IData histogramData = histograms.get(i);
                    histogramSinks.get(i).putData(histogramData);
                }
            }
        }
        if (histogramList != null) {
            for (Histogram h : histogramList) {
                h.reset();
            }
        }
    }

    public IDataInfo processDataInfo(IDataInfo incomingDataInfo) {
        if (incomingDataInfo.getLength() > 1) {
            throw new RuntimeException("AccumulatorAverageCollapsingLog can only handle single data");
        }
        dataInfo = new DataInfoFunction(incomingDataInfo.getLabel(), incomingDataInfo.getDimension(), this);
        reset();
        return dataInfo;
    }

    public int getIndependentArrayDimension() {
        return 1;
    }

    public DataDoubleArray getIndependentData(int i) {
        return nData;
    }

    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return nDataInfo;
    }

    public DataTag getIndependentTag() {
        return nTag;
    }

    public void enableRawLogStatistics() {
        lAvg = new AccumulatorAverageCollapsing(200);
        lData = new DataDouble();
        lDataInfo = new DataInfoDouble("raw log", Null.DIMENSION);
        lDataInfo.addTag(new DataTag());
        lAvg.putDataInfo(lDataInfo);
    }

    public AccumulatorAverage getRawLogAccumulator() {
        return lAvg;
    }
}
