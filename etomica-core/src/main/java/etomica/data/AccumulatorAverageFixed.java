/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.math.function.Function;
import etomica.math.function.IFunction;

/**
 * AccumulatorAverage that maintains a fixed block size.
 * <p>
 * This accumulator accepts any type of Data
 */
public class AccumulatorAverageFixed extends AccumulatorAverage {

    protected final IFunction negativeChop, sanityCheckBC;
    protected IData sum; //sum(value)
    protected IData sumBlockSquare;//sum(blockAvg^2)
    protected IData currentBlockSum;//block_sum(value)
    protected IData sumSquare;//sum(value^2)
    protected IData mostRecentBlock, correlationSum, firstBlock;
    protected IData work, work2;
    protected boolean doStrictBlockData = false;
    protected IDataSink blockDataSink;

    /**
     * Default constructor sets block size to 1000 and sets the
     * interval for pushing the output data (pushInterval) to 100.
     */
    public AccumulatorAverageFixed() {
        this(1000);
    }

    /**
     * Constructs with interval for pushing the output data set to 100.
     *
     * @param blockSize size of the block.
     */
    public AccumulatorAverageFixed(long blockSize) {
        super(blockSize);
        // this chops negative values to be 0, used in cases where value
        // is non-zero due to roundoff and happens to be negative
        negativeChop = new IFunction() {
            public double f(double x) {
                return !(x <= 0) ? x : 0;
            }
        };

        // this provides a sanity check for the block correlation, which must
        // be between -1 and 1.  It can take strange (nonsense) values,
        // especially when the value is consistent (stdev=0).
        sanityCheckBC = new IFunction() {
            public double f(double x) {
                return (Double.isNaN(x) || x <= -1 || x >= 1) ? 0 : x;
            }
        };
    }

    /**
     * @return value of DoStrictBlockData
     */
    public boolean getDoStrictBlockData() {
        return doStrictBlockData;
    }

    /**
     * Sets the accumulator to base statistics only on block data (where that
     * makes sense).  The primary effect here is that data does not go into
     * the average unless it is part of a complete block.  Because of how it
     * is computed, the standard deviation is unaffected by this flag.
     *
     * @param newDoStrictBlockData
     */
    public void setDoStrictBlockData(boolean newDoStrictBlockData) {
        doStrictBlockData = newDoStrictBlockData;
    }

    /**
     * @return the data sink where block data are passed.
     */
    public IDataSink getBlockDataSink() {
        return blockDataSink;
    }

    /**
     * Allows individual block data to be passed to a data sink. Each block is
     * pushed to the sink after it is complete (as determined by block size).
     *
     * @param blockDataSink the sink where the blocks are to be passed. For
     *                      multiple streams, this could be a DataFork.
     */
    public void setBlockDataSink(IDataSink blockDataSink) {
        this.blockDataSink = blockDataSink;
        if (dataInfo != null) {
            blockDataSink.putDataInfo(((DataInfoGroup) dataInfo).getSubDataInfo(MOST_RECENT.index));
        }
    }

    /**
     * @return null (any data is good data)
     */
    public DataPipe getDataCaster(IEtomicaDataInfo incomingDataInfo) {
        return null;
    }

    /**
     * Add the given values to the sums and block sums. If any of the given data
     * values is NaN, method returns with no effect on accumulation sums.
     *
     * @param data data to be added.
     * @return true if data are not NaN.
     */
    public boolean addData(IData data) {
        if (data.isNaN())
            return false;

        mostRecent.E(data);
        currentBlockSum.PE(data);
        work.E(data);
        work.TE(data);
        sumSquare.PE(work);
        if (--blockCountDown == 0) {//count down to zero to determine
            // completion of block
            doBlockSum();
            if (blockDataSink != null) {
                blockDataSink.putData(mostRecentBlock);
            }
        }
        return true;
    }

    public void putDataInfo(IEtomicaDataInfo inputDataInfo) {
        super.putDataInfo(inputDataInfo);
        if (blockDataSink != null) {
            blockDataSink.putDataInfo(inputDataInfo);
        }
    }

    /**
     * Performs the block sum after <tt>blockSize</tt> calls to addData.
     */
    protected void doBlockSum() {
        count++;
        sum.PE(currentBlockSum);
        blockCountDown = blockSize;
        currentBlockSum.DE(blockSize);//compute block average
        work.E(currentBlockSum);
        work.TE(currentBlockSum);
        sumBlockSquare.PE(work);
        if (count > 1) {
            work.E(currentBlockSum);
            work.TE(mostRecentBlock);
            correlationSum.PE(work);
        } else {
            firstBlock.E(currentBlockSum);
        }
        //reset blocks
        mostRecentBlock.E(currentBlockSum);
        currentBlockSum.E(0.0);
    }

    public IData getData() {
        if (sum == null)
            return null;
        if (count > 0) {
            // calculate block average (discarded if !doStrictBlockData)
            average.E(sum);
            average.DE(count * blockSize);
            work.E(average);
            work.TE(average);
        }
        if (count > 1) {
            // calculate other block properties (these require 2 or more blocks)

            error.E(sumBlockSquare);
            error.DE(count);
            error.ME(work);
            error.map(negativeChop);

            // error's intermediate value is useful for calculating block correlation
            blockCorrelation.E(average);
            blockCorrelation.TE(-2 * count);
            blockCorrelation.PE(firstBlock);
            blockCorrelation.PE(mostRecentBlock);
            blockCorrelation.TE(average);
            blockCorrelation.PE(correlationSum);
            blockCorrelation.DE(count - 1);
            blockCorrelation.PE(work);
            blockCorrelation.DE(error);
            blockCorrelation.map(sanityCheckBC);

            // ok, now finish up with error
            error.DE(count - 1);
            error.map(Function.Sqrt.INSTANCE);

            if (doIncludeACInError && count > 3) {
                work.E(1); // 1+c
                work.PE(blockCorrelation);
                work2.E(1);
                work2.ME(blockCorrelation); // 1-c
                work.DE(work2); // (1+c)/(1-c)
                work.map(Function.Sqrt.INSTANCE);
                error.TE(work);
            }
        } else {
            error.E(Double.NaN);
            blockCorrelation.E(Double.NaN);
        }

        long nTotalData = count * blockSize + (blockSize - blockCountDown);
        if (nTotalData > 0) {
            // now use *all* of the data
            if (!doStrictBlockData) {
                average.E(sum);
                average.PE(currentBlockSum);
                average.DE(nTotalData);
                work.E(average);
                work.TE(average);
            } else {
                work.E(sum);
                work.PE(currentBlockSum);
                work.DE(nTotalData);
                work.TE(average);
            }
            standardDeviation.E(sumSquare);
            standardDeviation.DE(nTotalData);
            standardDeviation.ME(work);
            standardDeviation.map(negativeChop);
            standardDeviation.map(Function.Sqrt.INSTANCE);
        } else {
            average.E(Double.NaN);
            standardDeviation.E(Double.NaN);
        }
        return dataGroup;
    }

    public void reset() {
        super.reset();
        if (sum == null) {
            return;
        }
        sum.E(0);
        sumBlockSquare.E(0);
        currentBlockSum.E(0);
        sumSquare.E(0);
        correlationSum.E(0);
        firstBlock.E(Double.NaN);
        mostRecentBlock.E(Double.NaN);
    }

    public IEtomicaDataInfo processDataInfo(IEtomicaDataInfo incomingDataInfo) {
        sum = incomingDataInfo.makeData();
        sumBlockSquare = incomingDataInfo.makeData();
        currentBlockSum = incomingDataInfo.makeData();
        sumSquare = incomingDataInfo.makeData();
        firstBlock = incomingDataInfo.makeData();
        correlationSum = incomingDataInfo.makeData();
        mostRecentBlock = incomingDataInfo.makeData();
        work = incomingDataInfo.makeData();
        work2 = incomingDataInfo.makeData();
        return super.processDataInfo(incomingDataInfo);
    }
}
