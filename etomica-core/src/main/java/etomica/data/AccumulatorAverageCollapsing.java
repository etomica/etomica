/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;


import java.util.Arrays;

/**
 * AccumulatorAverage that adjusts the block size during the simulation.  When
 * a certain number of blocks of data have been collected, the blocks are
 * collapsed such that each new block contains the data from two of the old
 * blocks. In other words, with each collapse, the number of blocks currently holding
 * data is cut in half and the size of each block is doubled.
 * <p>
 * This allows the accumulator to yield reasonable estimates of the
 * average and uncertainty during the beginning of the simulation and still
 * produce accurate estimates of the uncertainty for longer simulation runs.
 * <p>
 * This accumulator can only operate on Data with a single value.
 */
public class AccumulatorAverageCollapsing extends AccumulatorAverage {

    protected int maxBlocks;
    protected double[] blockAvgs;
    protected double currentBlockAvg, totalSumSquare;
    protected double totalBlockVarSum, totalVarSum;
    protected double correlationSum;
    protected double totalBlockAvg, totalAvg;
    protected long nominalBlockSize;

    /**
     * Default constructor sets the maximum number of blocks (maxBlocks) to a default value of 20,
     * the initial block size (blockSize) to 1 and the interval for pushing the output data (pushInterval)
     * to 100.
     */
    public AccumulatorAverageCollapsing() {
        this(20);
    }

    /**
     * Constructor with default pushInterval of 100 and initial blockSize of 1.
     *
     * @param maxBlocks the largest number of blocks. Adding new data beyond this triggers
     *                  collapsing of blocks
     */
    public AccumulatorAverageCollapsing(int maxBlocks) {
        this(maxBlocks, 1);
    }

    /**
     * Constructor with default pushInterval of 100.
     *
     * @param maxBlocks the largest number of blocks. Adding new data beyond this triggers
     *                  collapsing of blocks
     * @param blockSize the initial number of terms that contribute to each block. Doubles
     *                  with each collapse of the blocks.
     */
    public AccumulatorAverageCollapsing(int maxBlocks, int blockSize) {
        super(blockSize);
        blockAvgs = new double[0];
        setMaxBlocks(maxBlocks);
        setBlockSize(blockSize);
        setPushInterval(100);
    }

    /**
     * @return the current value of maxBlocks.
     */
    public int getMaxBlocks() {
        return maxBlocks;
    }

    /**
     * Changes the maximum number of blocks. If the new value is less than the
     * current number of blocks, collapsing is done until number of blocks is
     * less than or equal to the given value.
     *
     * @param newMaxBlocks new value for maxBlocks
     * @throws RuntimeException if newMaxBlocks is less than 4 or an odd number.
     */
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
        blockAvgs = Arrays.copyOf(blockAvgs, maxBlocks);
    }

    public IDataInfo processDataInfo(IDataInfo incomingDataInfo) {
        if (incomingDataInfo.getLength() > 1) {
            // Can only handle single data.
            throw new RuntimeException("AccumulatorAverageCollapsing can only handle single data");
        }
        return super.processDataInfo(incomingDataInfo);
    }

    /**
     * Add the given values to the sums and block sums. If any of the given data
     * values is NaN, method returns with no effect on accumulation sums.
     *
     * @param data Data to be added.
     */
    public boolean addData(IData data) {
        if (data.isNaN())
            return false;
        double value = data.getValue(0);
        currentBlockAvg += (value - currentBlockAvg) / (blockSize - blockCountDown + 1);
        mostRecent.E(value);
        double t = totalAvg;
        totalAvg += (value - totalAvg) / (count * blockSize + (blockSize - blockCountDown + 1));
        totalVarSum += (value - t) * (value - totalAvg);

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
        double oldTotalBlockAvg = totalBlockAvg;
        totalBlockAvg += (currentBlockAvg - totalBlockAvg) / (count + 1);
        if (count > 0) {
            correlationSum += blockAvgs[(int) count - 1] * currentBlockAvg;
        }
        blockAvgs[(int) count] = currentBlockAvg;
        totalBlockVarSum += (currentBlockAvg - oldTotalBlockAvg) * (currentBlockAvg - totalBlockAvg);
        currentBlockAvg = 0;

        count++;
        blockCountDown = blockSize;

        if (count == maxBlocks) {
            collapseBlocks();
        }

    }

    protected void collapseBlocks() {
        int intCount = (int) count;
        if (intCount % 2 == 1) {
            // if we have an odd number of blocks, the last block will get
            // dropped.  So add its contribution to the "current" block.
            currentBlockAvg += blockAvgs[intCount - 1] * blockSize / (2 * blockSize - blockCountDown);
        }
        totalBlockVarSum = 0;
        count /= 2;
        intCount = (int) count;
        correlationSum = 0;
        // the first half of the blocks contain all previous data
        totalBlockAvg = 0;
        for (int i = 0; i < intCount; i++) {
            double t = totalBlockAvg;
            blockAvgs[i] = (blockAvgs[2 * i] + blockAvgs[2 * i + 1]) / 2;
            totalBlockAvg += (blockAvgs[i] - totalBlockAvg) / (i + 1);
            totalBlockVarSum += (blockAvgs[i] - t) * (blockAvgs[i] - totalBlockAvg);
            if (i > 0) {
                correlationSum += blockAvgs[i - 1] * blockAvgs[i];
            }
        }
        // the last half are 0
        for (int i = intCount; i < maxBlocks; i++) {
            blockAvgs[i] = 0;
        }
        blockCountDown += blockSize;
        blockSize *= 2;
    }

    /**
     * @return a DataGroup with Data instances holding the statistics kept by
     * this accumulator (as described in general comments for this class).
     */
    public IData getData() {
        if (dataGroup == null) {
            return null;
        }

        if (count > 0) {
            average.E(totalBlockAvg);
        } else {
            average.E(Double.NaN);
        }
        if (count > 1) {
            // calculate block properties (these require 2 or more blocks)
            double err = totalBlockVarSum / count;
            error.E(Math.sqrt(err / (count - 1)));

            double bc = ((correlationSum - (totalBlockAvg * count - blockAvgs[0]) * (totalBlockAvg * count - blockAvgs[(int) count - 1]) / (count - 1)) / (count - 1)) / err;
            // sanity check
            bc = (Double.isNaN(bc) || bc <= -1 || bc >= 1) ? 0 : bc;
            blockCorrelation.E(bc);

            if (doIncludeACInError && count > 3) {
                double c = blockCorrelation.getValue(0);
                error.TE(Math.sqrt((1 + c) / (1 - c)));
            }
        } else {
            error.E(Double.NaN);
            blockCorrelation.E(Double.NaN);
        }

        long nTotalData = count * blockSize + (blockSize - blockCountDown);
        if (nTotalData > 1) {
            standardDeviation.E(Math.sqrt(totalVarSum / (count * blockSize + (blockSize - blockCountDown))));
        } else {
            standardDeviation.E(Double.NaN);
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
        currentBlockAvg = 0;
        totalSumSquare = 0;
        totalBlockAvg = 0;
        totalBlockVarSum = 0;
        totalVarSum = 0;
    }

    /**
     * Sets the size of the block used to group data for error analysis. Resets
     * statistics accumulated so far.  Default is 1.
     *
     * @param newBlockSize new block size.
     */
    public void setBlockSize(long newBlockSize) {
        nominalBlockSize = newBlockSize;
        super.setBlockSize(newBlockSize);
    }
}
