/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.util.random.IRandom;
import etomica.util.random.RandomMersenneTwister;
import etomica.util.random.RandomNumberGeneratorUnix;

import java.util.Arrays;

public class AccumulatorAverageCollapsingLogAB extends
        AccumulatorAverageCollapsingLog {

    public AccumulatorAverageCollapsingLogAB() {
        this(new RandomMersenneTwister(RandomNumberGeneratorUnix.getRandSeedArray()));
    }

    public AccumulatorAverageCollapsingLogAB(IRandom random) {
        super(random);
        setMaxSample(1);
    }

    public void setMaxSample(double newMaxSample) {
        maxSample = newMaxSample;
    }
    
    protected void resizeData(int newSize) {
        super.resizeData(newSize);
        labSums = Arrays.copyOf(labSums, newSize);
        labSums[newSize-1] = new double[nMoments];
    }


    /**
     * Add the given values to the sums and block sums. If any of the given data
     * values is NaN, method returns with no effect on accumulation sums.
     */
    public boolean addData(IData data) {
        if (data.isNaN())
            return false;
        double value = data.getValue(0);
        if (value == 0 || value == Double.POSITIVE_INFINITY) {
            throw new RuntimeException("got "+value);
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
                    double l = Math.log(rawData[i]);
                    doSums(lSums[0], l);
                    doSums(labSums[0], 0.5*(Math.log(0.5*(rawData[i]+maxSample))+l));
                    doSums(sums[0], rawData[i+1]);
                    l = Math.log(rawData[i+1]);
                    doSums(lSums[0], l);
                    doSums(labSums[0], 0.5*(Math.log(0.5*(rawData[i+1]+maxSample))+l));
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
        double l = Math.log(value);
        doSums(lSums[0], l);
        doSums(labSums[0], 0.5*(Math.log(0.5*(value+maxSample))+l));
        blockSums[1] += value;
        count++;
        int nDoubles = Long.numberOfTrailingZeros(count);
        for (int i=1; i<nDoubles+1; i++) {
            long n = 1L<<i;
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
                        l = Math.log(rawData[j]);
                        doSums(lSums[i], l);
                        doSums(labSums[i], 0.5*(Math.log((rawData[j]*n+maxSample)/(n+1))+l));
                        doSums(sums[i], rawData[j+1]);
                        l = Math.log(rawData[j+1]);
                        doSums(lSums[i], l);
                        doSums(labSums[i], 0.5*(Math.log((rawData[j+1]*n+maxSample)/(n+1))+l));
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
            l = Math.log(blockSums[i]);
            doSums(lSums[i], l);
            doSums(labSums[i], 0.5*(Math.log((blockSums[i]*n+maxSample)/(n+1))+l));
            blockSums[i+1] += blockSums[i];
            blockSums[i] = 0;
        }
        return true;
    }

    /**
     * Returns the average natural logs of the block averages, for all block
     * sizes.
     */
    public IData getAverageAntibiasedLogs() {
        double[] av = averages.getData();
        if (av.length == 0) {
            return averages;
        }
        if (initialSeed >= 0) {
            ((RandomMersenneTwister)random).setSeed(initialSeed);
        }
        for (int i=0; i<rawDataBlockSize; i++) {
            av[i] = labSums[i][0] / (count>>i);
        }
        double labSum = 0;
        long n = 1L<<(rawDataBlockSize-1);
        for (int j=0; j<nRawData; j++) {
            labSum += 0.5*(Math.log(rawData[j]) + Math.log((rawData[j]*n+maxSample)/(n+1)));
            rawData2[j] = rawData[j];
        }
        av[rawDataBlockSize] = labSum/nRawData;
        for (int i=rawDataBlockSize+1; i<rawDataBlockSize+(32-Integer.numberOfLeadingZeros(nRawData)); i++) {
            labSum = 0;
            n *= 2;
            if (withReplacement) {
                reblockData(1<<(i-rawDataBlockSize));
            }
            else {
                reblockData();
            }
            for (int j=0; j<nRawData; j++) {
                labSum += 0.5*(Math.log(rawData2[j])+Math.log((rawData2[j]*n+maxSample)/(n+1)));
            }
            av[i] = labSum/nRawData;
        }
        return averages;
    }
    
    /**
     * Returns the standard deviation of the log of the block averages, for all
     * block sizes.
     */
    public IData getStdevantibiasedLog() {
        if (initialSeed >= 0) {
            ((RandomMersenneTwister)random).setSeed(initialSeed);
        }
        double[] std = stdev.getData();
        if (std.length == 0) {
            return stdev;
        }
        for (int i=0; i<rawDataBlockSize; i++) {
            long iCount = count >> i;
            double av = labSums[i][0] / iCount;
            std[i] = Math.sqrt(labSums[i][1] / iCount - av*av);
        }
        double iSum = 0, iSum2 = 0;
        long n = 1L<<(rawDataBlockSize-1);
        for (int j=0; j<nRawData; j++) {
            double l = 0.5*(Math.log(rawData[j])+Math.log((rawData[j]*n+maxSample)/(n+1)));
            iSum += l;
            iSum2 += l*l;
            rawData2[j] = rawData[j];
        }
        std[rawDataBlockSize] = Math.sqrt(iSum2/nRawData - (iSum*iSum)/(nRawData*nRawData));
        for (int i=rawDataBlockSize+1; i<rawDataBlockSize+(32-Integer.numberOfLeadingZeros(nRawData)); i++) {
            n *= 2;
            iSum = 0;
            iSum2 = 0;
            if (withReplacement) {
                reblockData(1<<(i-rawDataBlockSize));
            }
            else {
                reblockData();
            }
            for (int j=0; j<nRawData; j++) {
                double l = 0.5*(Math.log(rawData2[j])+Math.log((rawData2[j]*n+maxSample)/(n+1)));
                iSum += l;
                iSum2 += l*l;
            }
            std[i] = Math.sqrt(iSum2/nRawData - (iSum*iSum)/(nRawData*nRawData));
        }
        return stdev;
    }


    
    /**
     * Resets all sums to zero. All statistics are cleared.
     */
    public void reset() {
        super.reset();
        labSums = new double[0][nMoments];
    }


    protected double[][] labSums;
    protected double maxSample;
}
