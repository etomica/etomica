package etomica.data;

import etomica.util.Function;

/**
 * AccumulatorAverage that maintains a fixed block size.
 * <p>
 * This accumulator accepts any type of Data
 */
public class AccumulatorAverageFixed extends AccumulatorAverage {

    /**
     * Default constructor sets block size to Default value, and sets the
     * interval for pushing the output data (pushInterval) to 100.
     */
    public AccumulatorAverageFixed() {
        this(1000);
    }
    
    public AccumulatorAverageFixed(int blockSize) {
        super(blockSize);
    }

    /**
     * Returns null (any data is good data)
     */
    public DataPipe getDataCaster(IDataInfo incomingDataInfo) {
        return null;
    }

    /**
     * Add the given values to the sums and block sums. If any of the given data
     * values is NaN, method returns with no effect on accumulation sums.
     */
    public void addData(Data data) {
        if (data.isNaN())
            return;

        mostRecent.E(data);
        blockSum.PE(data);
        work.E(data);
        work.TE(data);
        blockSumSq.PE(work);
        if (--blockCountDown == 0) {//count down to zero to determine
                                    // completion of block
            doBlockSum();
        }
    }
    
    /**
     * Performs the block sum after <tt>blockSize</tt> calls to addData.
     */
    protected void doBlockSum() {
        count++;
        blockCountDown = blockSize;
        blockSum.TE(1 / (double) blockSize);//compute block average
        //XXX should we divide blockSumSq by blockSize too?
        sum.PE(blockSum);
        work.E(blockSum);
        work.TE(blockSum);
        sumSquare.PE(work);
        if (!mostRecentBlock.isNaN()) {
            work.E(blockSum);
            work.TE(mostRecentBlock);
            correlationSum.PE(work);
        }
        else {
            firstBlock.E(blockSum);
        }
        sumSquareBlock.PE(blockSumSq);
        //reset blocks
        mostRecentBlock.E(blockSum);
        blockSum.E(0.0);
        blockSumSq.E(0.0);
    }

    /**
     * Returns a DataGroup with Data instances holding the statistics kept by
     * this accumulator (as described in general comments for this class).
     */
    public Data getData() {
        if (sum == null)
            return null;
        if (count > 0) {
            average.E(sum);
            average.TE(1 / (double) count);
            work.E(average);
            work.TE(average);
            error.E(sumSquare);
            error.TE(1 / (double) count);
            error.ME(work);

            // error's intermediate value is useful for calculating block correlation
            blockCorrelation.E(sum);
            blockCorrelation.TE(2.0);
            blockCorrelation.ME(firstBlock);
            blockCorrelation.ME(mostRecentBlock);
            blockCorrelation.TE(average);
            blockCorrelation.ME(correlationSum);
            blockCorrelation.TE(1.0/(1-count));
            blockCorrelation.PE(work);
            blockCorrelation.DE(error);
            
            // ok, now finish up with error
            error.TE(1 / (double) (count - 1));
            error.map(Function.Sqrt.INSTANCE);
            standardDeviation.E(sumSquareBlock);
            standardDeviation.TE(1.0 / (count * blockSize));
            standardDeviation.ME(work);
            standardDeviation.map(Function.Sqrt.INSTANCE);
        }
        else if (blockSize > blockCountDown) {
            average.E(blockSum);
            average.TE(1.0/(blockSize-blockCountDown));
            error.E(Double.NaN);
            blockCorrelation.E(Double.NaN);
            standardDeviation.E(Double.NaN);
        }
        return dataGroup;
    }

    /**
     * Resets all sums to zero. All statistics are cleared.
     */
    public void reset() {
        if (count == 0 && blockCountDown == blockSize) {
            //no data has been added yet, so nothing to reset
            return;
        }
        super.reset();
        sum.E(0);
        sumSquare.E(0);
        sumSquareBlock.E(0);
        blockSum.E(0);
        blockSumSq.E(0);
        correlationSum.E(0);
        firstBlock.E(Double.NaN);
        mostRecentBlock.E(Double.NaN);
    }

    /**
     * Sets the size of the block used to group data for error analysis. Has no
     * effect on statistics accumulated so far.  Default is 100.
     * 
     * @param blockSize
     *            new block size.
     */
    public void setBlockSize(int blockSize) {
        this.blockSize = blockSize;
        blockCountDown = blockSize;
    }

    /**
     * Returns the current value of the block size.
     */
    public int getBlockSize() {
        return blockSize;
    }

    /**
     * Prepares the accumulator for input data.  Discards any previous 
     * contributions to statistics.
     * 
     * @param incomingDataInfo
     *            the DataInfo instance for the data that will be given to
     *            addData
     */
    public IDataInfo processDataInfo(IDataInfo incomingDataInfo) {
        sum = incomingDataInfo.makeData();
        sumSquare = incomingDataInfo.makeData();
        sumSquareBlock = incomingDataInfo.makeData();
        blockSum = incomingDataInfo.makeData();
        blockSumSq = incomingDataInfo.makeData();
        firstBlock = incomingDataInfo.makeData();
        correlationSum = incomingDataInfo.makeData();
        mostRecentBlock = incomingDataInfo.makeData();
        work = incomingDataInfo.makeData();
        return super.processDataInfo(incomingDataInfo);
    }


    private static final long serialVersionUID = 1L;
    protected Data sum; //sum(blockSum/blkSize) = sum(blockAvg)
    protected Data sumSquare;//sum(blockAvg^2)
    protected Data blockSum;//block(value)
    protected Data blockSumSq;//block(value^2)
    protected Data sumSquareBlock;//sum(value^2)
    protected Data mostRecentBlock, correlationSum, firstBlock;
    protected Data work;
}
