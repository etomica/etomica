/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.util.Arrays;

/**
 * AccumulatorAverage that adjusts the block size during the simulation.  When
 * a certain number of blocks of data have been collected, the blocks are
 * collapsed such that each new block contains the data from two of the old
 * blocks.  This allows the accumulator to yield reasonable estimates of the
 * average and uncertainty during the beginning of the simulation and still
 * produce accurate estimates of the uncertainty for longer simulation runs.
 * <p>
 * This accumulator can only operate on Data with a single value.
 */
public class AccumulatorAverageCollapsing extends AccumulatorAverage {

    /**
     * Default constructor sets block size to Default value, and sets the
     * interval for pushing the output data (pushInterval) to 100.
     */
    public AccumulatorAverageCollapsing() {
        this(20);
    }
    
    public AccumulatorAverageCollapsing(int maxBlocks) {
        this(maxBlocks, 1);
    }
    
    public AccumulatorAverageCollapsing(int maxBlocks, int blockSize) {
        super(blockSize);
        blockSums = new double[0];
        setMaxBlocks(maxBlocks);
        setBlockSize(blockSize);
        setPushInterval(100);
    }

    public void setMaxBlocks(int newMaxBlocks) {
        // count is declared to be a long, but we know we can cast it to an int
        // because maxBlocks is an int
        if (newMaxBlocks < 4) {
            throw new RuntimeException("You need at least 3 blocks.  Get real!");
        }
        if (newMaxBlocks % 2 != 0) {
            throw new RuntimeException("I can only handle an even number of blocks");
        }
        while (count >= newMaxBlocks) {
            collapseBlocks();
        }
        maxBlocks = newMaxBlocks;
        // this would drop stuff at the end, but then we just collapsed the
        // blocks so they should all be 0.
        blockSums = Arrays.resizeArray(blockSums, maxBlocks);
    }
    
    public int getMaxBlocks() {
        return maxBlocks;
    }
    
    /**
     * Checks that incoming Data implements Data, and returns null if
     * this is so. Otherwise throws a ClassCastException, as there is no data
     * caster to Data.
     */
    public DataPipe getDataCaster(IEtomicaDataInfo incomingDataInfo) {
        if (incomingDataInfo.getLength() > 1) {
            throw new RuntimeException("AccumulatorAverageCollapsing can only handle single data");
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
        currentBlockSum += value;
        totalSumSquare += value * value;
        mostRecent.E(value);
        if (--blockCountDown == 0) {//count down to zero to determine
                                    // completion of block
            doBlockSum();
        }
        return true;
    }
    
    /**
     * Performs the block sum after <tt>blockSize</tt> calls to addData.
     */
    protected void doBlockSum() {
        
        totalBlockSum += currentBlockSum;
        currentBlockSum /= blockSize;
        if (count > 0) {
            correlationSum += blockSums[(int)count-1] * currentBlockSum;
        }
        blockSums[(int)count] = currentBlockSum;
        totalSumBlockSq += currentBlockSum * currentBlockSum;
        currentBlockSum = 0;
        
        count++;
        blockCountDown = blockSize;

        if (count == maxBlocks) {
            collapseBlocks();
        }

    }
    
    protected void collapseBlocks() {
        int intCount = (int)count;
        if (intCount % 2 == 1) {
            // if we have an odd number of blocks, the last block will get
            // dropped.  So add its contribution to the "current" block.
            currentBlockSum += blockSums[intCount-1];
            totalBlockSum -= blockSums[intCount-1];
        }
        totalSumBlockSq = 0;
        count /= 2;
        intCount = (int)count;
        correlationSum = 0;
        // the first half of the blocks contain all previous data
        for (int i=0; i<intCount; i++) {
            blockSums[i] = (blockSums[2*i] + blockSums[2*i+1]) / 2;
            totalSumBlockSq += blockSums[i] * blockSums[i];
            if (i>0) {
                correlationSum += blockSums[i-1]*blockSums[i];
            }
        }
        // the last half are 0
        for (int i=intCount; i<maxBlocks; i++) {
            blockSums[i] = 0;
        }
        blockCountDown += blockSize;
        blockSize *= 2;
    }

    /**
     * Returns a DataGroup with Data instances holding the statistics kept by
     * this accumulator (as described in general comments for this class).
     */
    public IData getData() {
    	if (dataGroup == null) {
            return null;
        }
    	
        if (count > 1) {
            // calculate block properties (these require 2 or more blocks)
            double blockAvg = totalBlockSum / (count*blockSize);
            double err = totalSumBlockSq / count - blockAvg * blockAvg;
            if (err < 0) {
                err = 0;
            }
            error.E(Math.sqrt(err / (count-1)));

            double bc = (((2 * totalBlockSum / blockSize - blockSums[0] - blockSums[(int)count-1]) * blockAvg - correlationSum) / 
                    (1-count) + blockAvg * blockAvg) / err;
            // sanity check
            bc = (Double.isNaN(bc) || bc <= -1 || bc >= 1) ? 0 : bc;
            blockCorrelation.E(bc);

            if (doIncludeACInError && count > 3) {
                double c = blockCorrelation.getValue(0);
                error.TE(Math.sqrt((1+c)/(1-c)));
            }
        }
        else {
            error.E(Double.NaN);
            blockCorrelation.E(Double.NaN);
        }

    	long nTotalData = count*blockSize + (blockSize-blockCountDown);
        if (nTotalData > 0) {
            double avg = (totalBlockSum+currentBlockSum) / nTotalData;
            average.E(avg);
            double variance = totalSumSquare / nTotalData - avg*avg;
            if (variance < 0) {
                variance = 0;
            }
            standardDeviation.E(Math.sqrt(variance));
        }
        return dataGroup;
    }

    /**
     * Resets all sums to zero. All statistics are cleared.
     */
    public void reset() {
        blockSize = nominalBlockSize;
        super.reset();

        correlationSum = 0;
        currentBlockSum = 0;
        totalSumSquare = 0;
        totalBlockSum = 0;
        totalSumBlockSq = 0;
    }

    /**
     * Sets the size of the block used to group data for error analysis. Resets
     * statistics accumulated so far.  Default is 1.
     * 
     * @param newBlockSize
     *            new block size.
     */
    public void setBlockSize(long newBlockSize) {
        nominalBlockSize = newBlockSize;
        super.setBlockSize(newBlockSize);
    }

    private static final long serialVersionUID = 1L;
    protected int maxBlocks;
    protected double[] blockSums;
    protected double currentBlockSum, totalSumSquare;
    protected double totalSumBlockSq;
    protected double correlationSum;
    protected double totalBlockSum;
    protected long nominalBlockSize;
}
