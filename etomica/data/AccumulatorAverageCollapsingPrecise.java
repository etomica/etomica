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
public class AccumulatorAverageCollapsingPrecise extends AccumulatorAveragePrecise {

    /**
     * Default constructor sets block size to Default value, and sets the
     * interval for pushing the output data (pushInterval) to 100.
     */
    public AccumulatorAverageCollapsingPrecise() {
        this(20);
    }
    
    public AccumulatorAverageCollapsingPrecise(int maxBlocks) {
        this(maxBlocks, 1);
    }
    
    public AccumulatorAverageCollapsingPrecise(int maxBlocks, int blockSize) {
        super(blockSize);
        blockAvgs = new double[0];
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
        blockAvgs = Arrays.resizeArray(blockAvgs, maxBlocks);
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
        count++;
        blockDataCount++;
        double value = data.getValue(0);
        currentBlockAvg += (value - currentBlockAvg) / blockDataCount;
        totalAvg += (value - totalAvg) / count;
        totalSquareAvg += (value * value - totalSquareAvg) / count;
        mostRecent.E(value);
        if (blockDataCount == blockSize) {//count down to zero to determine
                                    // completion of block
            doBlockSum();
        }
        return true;
    }
    
    /**
     * Performs the block sum after <tt>blockSize</tt> calls to addData.
     */
    protected void doBlockSum() {
        blockAvg = totalAvg;

        if (blockCount > 0) {
            correlationAvg += (blockAvgs[(int)blockCount-1] * currentBlockAvg - correlationAvg) / blockCount;
        }
        blockAvgs[(int)blockCount] = currentBlockAvg;
        totalSumBlockSq += currentBlockAvg * currentBlockAvg;
        currentBlockAvg = 0;
        
        blockCount++;
        blockDataCount = 0;

        if (blockCount == maxBlocks) {
            collapseBlocks();
        }
    }
    
    protected void collapseBlocks() {
        int intCount = (int)blockCount;
        if (intCount % 2 == 1) {
            // if we have an odd number of blocks, the last block will get
            // dropped.  So add its contribution to the "current" block.
            blockAvg += (blockAvg - blockAvgs[intCount-1]) / (blockCount-1);
            currentBlockAvg = (blockAvgs[intCount-1]*blockSize + currentBlockAvg*blockDataCount) / (blockSize+blockDataCount);
            blockDataCount += blockSize;
        }
        totalSumBlockSq = 0;
        blockCount /= 2;
        intCount = (int)blockCount;
        correlationAvg = 0;
        // the first half of the blocks contain all previous data
        for (int i=0; i<intCount; i++) {
            blockAvgs[i] = (blockAvgs[2*i] + blockAvgs[2*i+1]) / 2;
            totalSumBlockSq += blockAvgs[i] * blockAvgs[i];
            if (i>0) {
                correlationAvg += (blockAvgs[i-1]*blockAvgs[i] - correlationAvg) / i;
            }
        }
        // the last half are 0
        for (int i=intCount; i<maxBlocks; i++) {
            blockAvgs[i] = 0;
        }
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
    	
        if (blockCount > 1) {
            // calculate block properties (these require 2 or more blocks)
            double err = totalSumBlockSq / blockCount - blockAvg * blockAvg;
            if (err < 0) {
                err = 0;
            }
//            System.out.println("new block variance "+err);
            error.E(Math.sqrt(err / (blockCount-1)));

            blockCorrelation.E((correlationAvg - (2 * blockAvg * blockCount - blockAvgs[0] - blockAvgs[(int)blockCount-1]) / (blockCount-1) * blockAvg + blockAvg * blockAvg) / err);
        }
        else {
            error.E(Double.NaN);
            blockCorrelation.E(Double.NaN);
        }

        if (count > 0) {
            average.E(totalAvg);
            double variance = totalSquareAvg - totalAvg*totalAvg;
            if (variance < 0) {
                System.out.println("new negative variance "+variance);
                variance = 0;
            }
//            else {
//                System.out.println("new positive variance "+variance);
//            }
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

        blockDataCount = 0;
        correlationAvg = 0;
        currentBlockAvg = 0;
        totalSquareAvg = 0;
        totalAvg = 0;
        totalSumBlockSq = 0;
    }

    /**
     * Sets the size of the block used to group data for error analysis. Resets
     * statistics accumulated so far.  Default is 1.
     * 
     * @param blockSize
     *            new block size.
     */
    public void setBlockSize(long newBlockSize) {
        nominalBlockSize = newBlockSize;
        super.setBlockSize(newBlockSize);
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
    public void addDataSink(IDataSink newDataSink, StatType[] types) {
        int[] indexes = new int[types.length];
        for (int i = 0; i < types.length; i++) {
            indexes[i] = types[i].index;
        }
        DataGroupFilter filter = new DataGroupFilter(indexes);
        addDataSink(filter);
        filter.setDataSink(newDataSink);
    }

    private static final long serialVersionUID = 1L;
    protected int maxBlocks;
    protected double[] blockAvgs;
    protected double currentBlockAvg, totalSquareAvg;
    protected double totalSumBlockSq;
    protected double correlationAvg;
    protected double totalAvg, blockAvg;
    protected long nominalBlockSize;
    protected long blockDataCount;
}
