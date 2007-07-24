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
        blockSumSqs = new double[0];
        setMaxBlocks(maxBlocks);
        setBlockSize(blockSize);
        setPushInterval(100);
        
    }

    public void setMaxBlocks(int newMaxBlocks) {
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
        blockSumSqs = Arrays.resizeArray(blockSumSqs, maxBlocks);
    }
    
    public int getMaxBlocks() {
        return maxBlocks;
    }
    
    /**
     * Checks that incoming Data implements Data, and returns null if
     * this is so. Otherwise throws a ClassCastException, as there is no data
     * caster to Data.
     */
    public DataPipe getDataCaster(IDataInfo incomingDataInfo) {
        if (incomingDataInfo.getLength() > 1) {
            throw new RuntimeException("AccumulatorAverageCollapsing can only handle single data");
        }
        return null;
    }

    /**
     * Add the given values to the sums and block sums. If any of the given data
     * values is NaN, method returns with no effect on accumulation sums.
     */
    public void addData(Data data) {
        if (data.isNaN())
            return;
        double value = data.getValue(0);
        currentBlockSum += value;
        currentBlockSumSq += value * value;
        mostRecent.E(value);
        if (--blockCountDown == 0) {//count down to zero to determine
                                    // completion of block
            doBlockSum();
        }
    }
    
    /**
     * Performs the block sum after <tt>blockSize</tt> calls to addData.
     */
    protected void doBlockSum() {
        
        if (count > 0) {
            correlationSum *= blockSums[count-1] * currentBlockSum;
        }
        currentBlockSum /= blockSize;
        blockSums[count] = currentBlockSum;
        blockSumSqs[count] = currentBlockSumSq;
        totalSumBlockSq += currentBlockSum * currentBlockSum;
        totalSum += currentBlockSum;
        totalSumSq += currentBlockSumSq;
        currentBlockSum = 0;
        currentBlockSumSq = 0;
        
        count++;
        blockCountDown = blockSize;

        if (count == maxBlocks) {
            collapseBlocks();
        }

    }
    
    protected void collapseBlocks() {
        if (count % 2 == 1) {
            // if we have an odd number of blocks, the last block will get
            // dropped.  So add its contribution to the "current" block.
            currentBlockSum += blockSums[count-1];
            currentBlockSumSq += blockSumSqs[count-1];
            // we can just fix totalSum here rather than recalculating it
            totalSum -= blockSums[count-1];
            totalSumSq -= blockSumSqs[count-1];
        }
        totalSum /= 2;
        totalSumBlockSq = 0;
        count /= 2;
        // the first half of the blocks contain all previous data
        for (int i=0; i<count; i++) {
            blockSums[i] = (blockSums[2*i] + blockSums[2*i+1]) / 2;
            totalSumBlockSq += blockSums[i] * blockSums[i];
            blockSumSqs[i] = blockSumSqs[2*i] + blockSumSqs[2*i+1];
        }
        // the last half are 0
        for (int i=count; i<maxBlocks; i++) {
            blockSums[i] = 0;
            blockSumSqs[i] = 0;
        }
        blockCountDown += blockSize;
        blockSize *= 2;
    }

    /**
     * Returns a DataGroup with Data instances holding the statistics kept by
     * this accumulator (as described in general comments for this class).
     */
    public Data getData() {
    	if (dataGroup == null) {
            return null;
        }
        else if (count > 0) {
            double avg = totalSum / count;
            average.E(avg);
            double err = totalSumBlockSq / count - avg * avg;

            // error's intermediate value is useful for calculating block correlation
            blockCorrelation.E((((2 * totalSum - blockSums[0]) * avg - correlationSum) / (1-count) + avg * avg) / err);
            
            // ok, now finish up with error
            error.E(Math.sqrt(err / (count-1)));
            double stdev = Math.sqrt(totalSumSq / (count*blockSize) - avg*avg);
            standardDeviation.E(stdev);
        }
        else if (blockSize > blockCountDown) {
            average.E(currentBlockSum/(blockSize-blockCountDown));
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
        count = 0;
        blockSize = nominalBlockSize;
        blockCountDown = blockSize;
        correlationSum = 0;
        currentBlockSum = 0;
        currentBlockSumSq = 0;
        totalSum = 0;
        totalSumBlockSq = 0;
        totalSumSq = 0;

        error.E(Double.NaN);
        mostRecent.E(Double.NaN);
        average.E(Double.NaN);
        standardDeviation.E(Double.NaN);
        blockCorrelation.E(Double.NaN);
    }

    /**
     * Sets the size of the block used to group data for error analysis. Has no
     * effect on statistics accumulated so far.  Default is 100.
     * 
     * @param blockSize
     *            new block size.
     */
    public void setBlockSize(int newBlockSize) {
        nominalBlockSize = newBlockSize;
        blockSize = newBlockSize;
        blockCountDown = newBlockSize;
    }

    /**
     * Returns the current value of the block size.
     */
    public int getBlockSize() {
        return blockSize;
    }

    /**
     * Adds a new DataSink that will receive a specific subset of the statistics
     * generated by this accumulator. The DataSink will be given a DataGroup
     * that has only the specified statistics in it. The output of this
     * accumulator is passed through a DataGroupFilter that removes the unwanted
     * Data.
     * 
     * @param newDataSink
     *            the new DataSink
     * @param types
     *            array indicating the statistics to be included in the
     *            DataGroup sent to the sink.
     */
    public void addDataSink(DataSink newDataSink, StatType[] types) {
        int[] indexes = new int[types.length];
        for (int i = 0; i < types.length; i++) {
            indexes[i] = types[i].index;
        }
        DataGroupFilter filter = new DataGroupFilter(indexes);
        addDataSink(filter);
        filter.setDataSink(newDataSink);
    }

    /**
     * Returns the number of times addData has been called since construction or
     * the last call to reset.
     */
    public int getCount() {
        return count;
    }

    private static final long serialVersionUID = 1L;
    protected int maxBlocks;
    protected double[] blockSums, blockSumSqs;
    protected double currentBlockSum, currentBlockSumSq;
    protected double totalSumBlockSq;
    protected double correlationSum;
    protected double totalSum;
    protected double totalSumSq;
    protected int nominalBlockSize;
}
