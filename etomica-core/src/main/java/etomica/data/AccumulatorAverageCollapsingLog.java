/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import java.util.ArrayList;

import etomica.util.random.IRandom;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.units.Null;
import etomica.units.Quantity;
import etomica.util.Arrays;
import etomica.util.random.RandomMersenneTwister;
import etomica.util.random.RandomNumberGenerator;
import etomica.util.random.RandomNumberGeneratorUnix;


/**
 * AccumulatorAverage that keeps track of block averages and standard
 * deviations corresponding to all block sizes from 1 to the n where n is the
 * largest power of 2 less than the current number of samples.  In addition to
 * block averages and standard deviations of the sampled quantity, the natural
 * log of the block averages and standard deviation of the natural log is also
 * calculated (this is what makes this class actually useful).  getData returns
 * the average natural log values for all block sizes as a DataFunction.
 * <p>
 * This class also uses a method to improve the block statistics for large.
 * All block averages are kept for a particular block size (n1 blocks of size
 * b1).  Statistics for that block are computed from that data.  Block
 * averages for blocks of the next larger size (2*n1) are computed by randomly
 * selecting n1/2 pairs of block values from the b1 set (without including any
 * block value more than once).  After using all block values, n1/2 additional
 * pairs are taken from the full b1 set, such that n1 blocks of size b1*2 are
 * obtained, instead of n1/2 blocks that would normally exist.  Statistics for
 * block size b1*2 are computed from these n1 blocks.  This process is
 * repeated for each block size such that the largest block size (which would
 * normally have 1 block) has n1 blocks.
 * <p>
 * This accumulator can only operate on Data with a single value.
 */
public class AccumulatorAverageCollapsingLog extends DataAccumulator implements DataSourceIndependent {

    public AccumulatorAverageCollapsingLog() {
        this(2);
    }
    
    /**
     * With this constructor, the accumulator will use an internal RNG with a
     * set seed such that each resampling of the data will be identical (the
     * same set of samples will go into the averages and the standard
     * deviations, etc).
     * @param nMoments the number of moments handled by the accumulator. 
     *           should be 2 (standard deviation) or 3 (stdev + skew).
     */
    public AccumulatorAverageCollapsingLog(int nMoments) {
        this(new RandomMersenneTwister(RandomNumberGeneratorUnix.getRandSeedArray()), nMoments);
        initialSeed = random.nextInt(1<<30);
    }
    
    public AccumulatorAverageCollapsingLog(IRandom random) {
        this(random, 2);
    }
    
    /**
     * @param nMoments the number of moments handled by the accumulator. 
     *           should be 2 (standard deviation) or 3 (stdev + skew).
     */
    public AccumulatorAverageCollapsingLog(IRandom random, int nMoments) {
        if (nMoments < 2 || nMoments > 3) {
            throw new IllegalArgumentException("nMoments should be 2 or 3");
        }
        this.nMoments = nMoments;
        this.random = random;
        initialSeed = -1;
        setPushInterval(100);
        setNumRawDataDoubles(14);
        nTag = new DataTag();
    }

    /**
     * Sets the number of times the random block generation method is used.
     * The number of blocks used will be 2^numRawDataDoubles.
     */
    public void setNumRawDataDoubles(int newNumRawDataDoubles) {
        nRawDataDoubles = newNumRawDataDoubles;
        intArrayList = new ArrayList<Integer>(1<<nRawDataDoubles);
        reset();
    }

    /**
     * Checks that incoming Data implements Data, and returns null if
     * this is so. Otherwise throws a ClassCastException, as there is no data
     * caster to Data.
     */
    public DataPipe getDataCaster(IEtomicaDataInfo incomingDataInfo) {
        if (incomingDataInfo.getLength() > 1) {
            throw new RuntimeException("AccumulatorAverageCollapsingLog can only handle single data");
        }
        return null;
    }

    public IData processData(IData inputData) {
        if(!active) return null;
        int oldSize = sums.length;
        addData(inputData);
        if (sums.length > oldSize) {
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
            throw new RuntimeException("got "+value);
        }
        if (value == 0) {
            value = Double.MIN_VALUE;
        }
        if (lAvg != null) {
            lData.x = Math.log(value);
            lAvg.putData(lData);
        }
        if (rawDataBlockSize == 0) {
            // special case early part where we keep every sample as raw data
            rawData[nRawData] = value;
            nRawData++;
            count++;
            int nDoubles = 63-Long.numberOfLeadingZeros(count);
            if (nDoubles >= averages.getLength()) {
                resizeData(nDoubles+1);
            }
            if (nRawData == 1<<nRawDataDoubles) {
                // we filled our raw data.  do summation for blockSize=0,
                // prepare for blockSize=1
                for (int i=0; i<nRawData; i+=2) {
                    doSums(sums[0], rawData[i]);
                    doSums(lSums[0], Math.log(rawData[i]));
                    doSums(sums[0], rawData[i+1]);
                    doSums(lSums[0], Math.log(rawData[i+1]));
                    rawData[i/2] = 0.5*(rawData[i]+rawData[i+1]);
                }
                for (int i=nRawData/2; i<nRawData; i++) {
                    rawData[i] = 0;
                }
                rawDataBlockSize++;
                nRawData /= 2;
            }
            return true;
        }
        // process sample as blockSize=0 data
        doSums(sums[0], value);
        doSums(lSums[0], Math.log(value));
        blockSums[1] += value;
        count++;
        int nDoubles = Long.numberOfTrailingZeros(count);
        for (int i=1; i<nDoubles+1; i++) {
            blockSums[i] *= 0.5;
            if (rawDataBlockSize == i) {
                // we consider this block average to be raw data
                rawData[nRawData] = blockSums[i];
                blockSums[i] = 0;
                nRawData++;
                if (nRawData == 1<<nRawDataDoubles) {
                    // we filled our raw data
                    resizeData(nDoubles+1);
                    for (int j=0; j<nRawData; j+=2) {
                        doSums(sums[i], rawData[j]);
                        doSums(lSums[i], Math.log(rawData[j]));
                        doSums(sums[i], rawData[j+1]);
                        doSums(lSums[i], Math.log(rawData[j+1]));
                        rawData[j/2] = 0.5*(rawData[j]+rawData[j+1]);
                    }
                    for (int j=nRawData/2; j<nRawData; j++) {
                        rawData[j] = 0;
                    }
                    rawDataBlockSize++;
                    nRawData /= 2;
                }
                return true;
            }
            // process block average normally
            doSums(sums[i], blockSums[i]);
            doSums(lSums[i], Math.log(blockSums[i]));
            blockSums[i+1] += blockSums[i];
            blockSums[i] = 0;
        }
        return true;
    }

    public long getCount() {
        return count;
    }

    protected void resizeData(int newSize) {
        blockSums = Arrays.resizeArray(blockSums, newSize+1);
        sums = (double[][])Arrays.resizeArray(sums, newSize);
        sums[newSize-1] = new double[nMoments];
        lSums = (double[][])Arrays.resizeArray(lSums, newSize);
        lSums[newSize-1] = new double[nMoments];
        averages = new DataFunction(new int[]{newSize});
        stdev = new DataFunction(new int[]{newSize});
        if (nMoments > 2) {
            skew = new DataFunction(new int[]{newSize});
        }
        dataInfo = new DataInfoFunction(dataInfo.getLabel(), dataInfo.getDimension(), this);
        nDataInfo = new DataInfoDoubleArray("block size", Quantity.DIMENSION, new int[]{newSize});
        nDataInfo.addTag(nTag);
        nData = new DataDoubleArray(newSize);
        double[] n = nData.getData();
        for (int i=0; i<n.length; i++) {
            n[i] = 1L<<i;
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

    /**
     * Returns the averages for all block sizes
     */
    public IData getAverages() {
        double[] av = averages.getData();
        for (int i=0; i<rawDataBlockSize; i++) {
            av[i] = sums[i][0] / (count>>i);
        }
        int n = averages.getLength();
        double iSum = 0;
        for (int j=0; j<nRawData; j++) {
            iSum += rawData[j];
        }
        for (int i=rawDataBlockSize; i<n; i++) {
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
            ((RandomMersenneTwister)random).setSeed(initialSeed);
        }
        for (int i=0; i<rawDataBlockSize; i++) {
            av[i] = lSums[i][0] / (count>>i);
        }
        double lSum = 0;
        for (int j=0; j<nRawData; j++) {
            lSum += Math.log(rawData[j]);
            rawData2[j] = rawData[j];
        }
        av[rawDataBlockSize] = lSum/nRawData;
        for (int i=rawDataBlockSize+1; i<rawDataBlockSize+(32-Integer.numberOfLeadingZeros(nRawData)); i++) {
            lSum = 0;
            if (withReplacement) {
                reblockData(1<<(i-rawDataBlockSize));
            }
            else {
                reblockData();
            }
            for (int j=0; j<nRawData; j++) {
                lSum += Math.log(rawData2[j]);
            }
            av[i] = lSum/nRawData;
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
            for (int j=0; j<nRawData; j++) {
                ensureHappyIntArrayList();
                int idxidx1 = random.nextInt(intArrayList.size());
                int idx1 = intArrayList.get(idxidx1);
                intArrayList.remove(idxidx1);
                ensureHappyIntArrayList();
                int idxidx2 = random.nextInt(intArrayList.size());
                int idx2 = intArrayList.get(idxidx2);
                intArrayList.remove(idxidx2);
                double sum12 = (rawData2[idx1] + rawData2[idx2])/2;
                rawData3[j] = sum12;
            }
        }
        else {
            // with replacement
            ensureHappyIntArrayList();
            for (int j=0; j<nRawData; j++) {
                int idx1 = random.nextInt(nRawData);
                int idx2 = random.nextInt(nRawData);
                rawData3[j] = 0.5*(rawData2[idx1] + rawData2[idx2]);
            }
        }
        for (int j=0; j<nRawData; j++) {
            rawData2[j] = rawData3[j];
        }
    }
    
    /**
     * Reblock data from rawData into rawData2, taking blocks of size n
     */
    protected void reblockData(int n) {
        // with replacement
        for (int j=0; j<nRawData; j++) {
            double innerSum = 0;
            for (int i=0; i<n; i++) {
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
            for (int j=0; j<nRawData; j++) {
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
            ((RandomMersenneTwister)random).setSeed(initialSeed);
        }
        double[] std = stdev.getData();
        if (std.length == 0) {
            return stdev;
        }
        for (int i=0; i<rawDataBlockSize; i++) {
            long iCount = count >> i;
            double av = sums[i][0] / iCount;
            std[i] = Math.sqrt(sums[i][1] / iCount - av*av);
            if (Double.isNaN(std[i])) {
                std[i] = 0;
            }
        }
        double iSum = 0, iSum2 = 0;
        for (int j=0; j<nRawData; j++) {
            iSum += rawData[j];
            iSum2 += rawData[j]*rawData[j];
            rawData2[j] = rawData[j];
        }
        iSum /= nRawData;
        iSum2 /= nRawData;
        double diff = iSum2 - iSum*iSum;
        if (diff < 0) {
            std[rawDataBlockSize] = 0;
        }
        else {
            std[rawDataBlockSize] = Math.sqrt(diff);
        }
        for (int i=rawDataBlockSize+1; i<rawDataBlockSize+(32-Integer.numberOfLeadingZeros(nRawData)); i++) {
            iSum = 0;
            iSum2 = 0;
            if (withReplacement) {
                reblockData(1<<(i-rawDataBlockSize));
            }
            else {
                reblockData();
            }
            for (int j=0; j<nRawData; j++) {
                iSum += rawData2[j];
                iSum2 += rawData2[j]*rawData2[j];
            }
            iSum /= nRawData;
            iSum2 /= nRawData;
            diff = iSum2 - iSum*iSum;
            if (diff < 0) {
                std[rawDataBlockSize] = 0;
            }
            else {
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
            ((RandomMersenneTwister)random).setSeed(initialSeed);
        }
        double[] std = stdev.getData();
        if (std.length == 0) {
            return stdev;
        }
        for (int i=0; i<rawDataBlockSize; i++) {
            long iCount = count >> i;
            double av = lSums[i][0] / iCount;
            std[i] = Math.sqrt(lSums[i][1] / iCount - av*av);
        }
        double iSum = 0, iSum2 = 0;
        for (int j=0; j<nRawData; j++) {
            double l = Math.log(rawData[j]);
            iSum += l;
            iSum2 += l*l;
            rawData2[j] = rawData[j];
        }
        std[rawDataBlockSize] = Math.sqrt(iSum2/nRawData - (iSum*iSum)/(nRawData*nRawData));
        for (int i=rawDataBlockSize+1; i<rawDataBlockSize+(32-Integer.numberOfLeadingZeros(nRawData)); i++) {
            iSum = 0;
            iSum2 = 0;
            if (withReplacement) {
                reblockData(1<<(i-rawDataBlockSize));
            }
            else {
                reblockData();
            }
            for (int j=0; j<nRawData; j++) {
                double l = Math.log(rawData2[j]);
                iSum += l;
                iSum2 += l*l;
            }
            std[i] = Math.sqrt(iSum2/nRawData - (iSum*iSum)/(nRawData*nRawData));
        }
        return stdev;
    }

    /**
     * Returns the standard deviation of the log of the block averages, for all
     * block sizes.
     */
    public IData getSkewLog() {
        if (initialSeed >= 0) {
            ((RandomMersenneTwister)random).setSeed(initialSeed);
        }
        if (nMoments < 3) {
            throw new RuntimeException("skew not handled");
        }
        double[] std = skew.getData();
        for (int i=0; i<rawDataBlockSize; i++) {
            long iCount = count >> i;
            double av = lSums[i][0] / iCount;
            double av2 = lSums[i][1] / iCount;
            double std3 = Math.pow(av2 - av*av, 1.5);
            std[i] = (lSums[i][2]/iCount - 3*av*av2 + 2*av*av*av)/std3;
        }
        double iSum = 0, iSum2 = 0, iSum3 = 0;
        for (int j=0; j<nRawData; j++) {
            double l = Math.log(rawData[j]);
            iSum += l;
            iSum2 += l*l;
            iSum3 += l*l*l;
            rawData2[j] = rawData[j];
        }
        double av = iSum / nRawData;
        double av2 = iSum2 / nRawData;
        double std3 = Math.pow(av2 - av*av, 1.5);
        std[rawDataBlockSize] = (iSum3/nRawData - 3*av*av2 + 2*av*av*av) / std3;
        int iMax = rawDataBlockSize+(32-Integer.numberOfLeadingZeros(nRawData));
        for (int i=rawDataBlockSize+1; i<iMax; i++) {
            iSum = 0;
            iSum2 = 0;
            iSum3 = 0;
            if (withReplacement) {
                reblockData(1<<(i-rawDataBlockSize));
            }
            else {
                reblockData();
            }
            for (int j=0; j<nRawData; j++) {
                double l = Math.log(rawData2[j]);
                iSum += l;
                iSum2 += l*l;
                iSum3 += l*l*l;
            }
            av = iSum / nRawData;
            av2 = iSum2 / nRawData;
            std3 = Math.pow(av2 - av*av, 1.5);
            std[i] = (iSum3/nRawData - 3*av*av2 + 2*av*av*av) / std3;
        }
        return skew;
    }

    /**
     * Performs the block sum after
     */
    protected void doSums(double[] subSums, double v) {
        double vpj = 1;
        for (int i=0; i<nMoments; i++) {
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

    /**
     * Returns the correlation between consecutive values of the log of the "raw"
     * data held by the accumulator.  Raw data may be individual samples, or
     * block averages, depending on how much data has been collected.  The
     * method returning the standard deviation of the log is based on the idea
     * that the "raw" data is uncorrelated.  If the correlation is small, the
     * standard deviation can be adjusted.  If large, then the standard
     * deviation is bogus (too small) and you must either collect more data
     * (so that the raw data are block averages with more samples) or decrease
     * the # of raw doubles (setNumRawDataDoubles) so that the accumulator blocks
     * data more aggressively.
     */
    public double getRawLogDataCorrelation() {
        if (nRawData < 3) return Double.NaN;
        double lSum = 0, lSum2 = 0, cSum = 0;
        for (int j=0; j<nRawData; j++) {
            double l = Math.log(rawData[j]);
            lSum += l;
            lSum2 += l*l;
            if (j>0) {
                cSum += l*Math.log(rawData[j-1]);
            }
        }
        double avg2 = (lSum/nRawData);
        avg2 *= avg2;
        double var = (lSum2/nRawData - avg2);
        if (var <= 0) return 0;
        
        double blockCorrelation = (((Math.log(rawData[0]) -2*lSum + Math.log(rawData[nRawData-1])) * (lSum/nRawData) + cSum)/(nRawData-1) + avg2)/var;
        blockCorrelation = (Double.isNaN(blockCorrelation) || blockCorrelation <= -1 || blockCorrelation >= 1) ? 0 : blockCorrelation;
        return blockCorrelation;
    }
    
    public double getSampleVariance() {
        if (rawDataBlockSize > 0) {
            return (count*sums[0][1] - sums[0][0]*sums[0][0])/(count*(count-1));
        }
        double sums0 = 0, sums1 = 0;
        for (int i=0; i<nRawData; i++) {
            sums0 += rawData[i];
            sums1 += rawData[i]*rawData[i];
        }
        return (count*sums1 - sums0*sums0)/(count*(count-1));
    }
    
    public double getSampleLogVariance() {
        if (rawDataBlockSize > 0) {
            return (count*lSums[0][1] - lSums[0][0]*lSums[0][0])/(count*(count-1));
        }
        double sums0 = 0, sums1 = 0;
        for (int i=0; i<nRawData; i++) {
            double l = Math.log(rawData[i]);
            sums0 += l;
            sums1 += l*l;
        }
        return (count*sums1 - sums0*sums0)/(count*(count-1));
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
        rawData = new double[1<<nRawDataDoubles];
        rawData2 = new double[1<<nRawDataDoubles];
        rawData3 = new double[1<<nRawDataDoubles];
        rawDataBlockSize = 0;
        nRawData = 0;
    }

    public IEtomicaDataInfo processDataInfo(IEtomicaDataInfo incomingDataInfo) {
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

    protected long count;
    protected double[] blockSums;
    protected double[][] sums, lSums;
    protected DataFunction averages, stdev, skew;
    protected DataInfoDoubleArray nDataInfo;
    protected DataDoubleArray nData;
    protected final DataTag nTag;
    protected double[] rawData, rawData2, rawData3;
    protected int nRawData, rawDataBlockSize;
    protected int nRawDataDoubles;
    protected ArrayList<Integer> intArrayList;
    protected final IRandom random;
    protected final int nMoments;
    protected int initialSeed;
    protected final boolean withReplacement = true;
    protected AccumulatorAverageCollapsing lAvg;
    protected DataDouble lData;
    protected DataInfoDouble lDataInfo;
    
    public static void main(String[] args) {
        AccumulatorAverageCollapsingLog ac = new AccumulatorAverageCollapsingLog(new RandomNumberGenerator());
        ac.putDataInfo(new DataInfoDouble("foo", Null.DIMENSION));
        DataDouble d = new DataDouble();
        RandomNumberGenerator rng = new RandomNumberGenerator();
        for (long i=0; i<1L<<25; i++) {
            d.x = Math.exp(rng.nextGaussian());
            ac.putData(d);
        }
        System.out.println(ac.getAverages().getValue(0)+" "+ac.getAverages().getValue(1)+" "+ac.getAverages().getValue(2)+" "+ac.getAverages().getValue(3));
        System.out.println(ac.getStdev().getValue(0)+" "+ac.getStdev().getValue(1)+" "+ac.getStdev().getValue(2)+" "+ac.getStdev().getValue(3));
        System.out.println(ac.getAverageLogs().getValue(0)+" "+ac.getAverageLogs().getValue(1)+" "+ac.getAverageLogs().getValue(2)+" "+ac.getAverageLogs().getValue(3));
        System.out.println(ac.getStdevLog().getValue(0)+" "+ac.getStdevLog().getValue(1)+" "+ac.getStdevLog().getValue(2)+" "+ac.getStdevLog().getValue(3));
    }
}
