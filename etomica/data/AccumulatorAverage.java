/*
 * History
 * Created on Jul 26, 2004 by kofke
 */
package etomica.data;

import etomica.data.types.DataArithmetic;
import etomica.data.types.DataGroup;
import etomica.data.types.DataArithmetic.DataInfoArithmetic;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.simulation.Simulation;
import etomica.util.EnumeratedType;
import etomica.util.Function;

/**
 * Accumulator that keeps statistics for averaging and error analysis. The
 * statisics are incremented with every call to addData. Statistics kept are
 * <ul>
 * <li>the value last given to addData
 * <li>the mean, or average
 * <li>the statistical error (67% confidence limits), obtained by standard
 * error analysis of block averages
 * <li>the standard deviation
 * <li>the most recent block average
 * </ul>
 * Each statistic is recorded by a Data instance of the same type as the
 * incoming data. The output Data is a DataGroup formed from these Data
 * instances (in the order given above).
 * <p>
 * A <i>block </i> is a subset of the input data, formed by grouping (for
 * example) 100 or 1000 successive contributions to addData. Averages for each
 * block are considered to be independent of other block averages, and the
 * confidence limits for the overall average is obtained as the standard error
 * of the mean of these block averages.
 * <p>
 * Incoming Data must implement DataArithmetic.
 */
public class AccumulatorAverage extends DataAccumulator {

    /**
     * Default constructor sets block size to Default value, and sets the
     * interval for pushing the output data (pushInterval) to 100.
     */
    public AccumulatorAverage(Simulation sim) {
        this(sim.getDefaults().blockSize);
    }
    
    public AccumulatorAverage(int blockSize) {
        super();
        setBlockSize(blockSize);
        setPushInterval(100);
        mostRecentTag = new DataTag();
        averageTag = new DataTag();
        errorTag = new DataTag();
        standardDeviationTag = new DataTag();
        mostRecentBlockTag = new DataTag();
        blockCorrelationTag = new DataTag();
    }

    public Object getTag(StatType statType) {
        if (statType == StatType.MOST_RECENT) {
            return mostRecentTag;
        }
        if (statType == StatType.AVERAGE) {
            return averageTag;
        }
        if (statType == StatType.ERROR) {
            return errorTag;
        }
        if (statType == StatType.STANDARD_DEVIATION) {
            return standardDeviationTag;
        }
        if (statType == StatType.MOST_RECENT_BLOCK) {
            return mostRecentBlockTag;
        }
        if (statType == StatType.BLOCK_CORRELATION) {
            return blockCorrelationTag;
        }
        return null;
    }

    /**
     * Checks that incoming Data implements DataArithmetic, and returns null if
     * this is so. Otherwise throws a ClassCastException, as there is no data
     * caster to DataArithmetic.
     */
    public DataProcessor getDataCaster(DataInfo incomingDataInfo) {
        if (incomingDataInfo instanceof DataInfoArithmetic) {
            return null;
        }
        throw new ClassCastException(
                "Data type cannot be handled by AccumulatorAverage");
    }

    /**
     * Add the given values to the sums and block sums. If any of the given data
     * values is NaN, method returns with no effect on accumulation sums.
     */
    public void addData(Data data) {
        DataArithmetic value = (DataArithmetic) data;
        if (value.isNaN())
            return;

        mostRecent.E(data);
        blockSum.PE(value);
        work.E(data);
        work.TE(value);
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
        return dataGroup;
    }

    /**
     * Resets all sums to zero. All statistics are cleared.
     */
    public void reset() {
        if (sum == null) {
            //no data has been added yet, so nothing to reset
            return;
        }
        count = 0;
        sum.E(0);
        sumSquare.E(0);
        sumSquareBlock.E(0);
        blockSum.E(0);
        blockSumSq.E(0);
        error.E(Double.NaN);
        mostRecent.E(Double.NaN);
        mostRecentBlock.E(Double.NaN);
        average.E(Double.NaN);
        standardDeviation.E(Double.NaN);
        blockCountDown = blockSize;
        blockCorrelation.E(Double.NaN);
        correlationSum.E(0);
        firstBlock.E(0);
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
    public DataInfo processDataInfo(DataInfo incomingDataInfo) {
        sum = (DataArithmetic)incomingDataInfo.makeData();
        sumSquare = (DataArithmetic)incomingDataInfo.makeData();
        sumSquareBlock = (DataArithmetic)incomingDataInfo.makeData();
        standardDeviation = (DataArithmetic)incomingDataInfo.makeData();
        average = (DataArithmetic)incomingDataInfo.makeData();
        error = (DataArithmetic)incomingDataInfo.makeData();
        blockSum = (DataArithmetic)incomingDataInfo.makeData();
        blockSumSq = (DataArithmetic)incomingDataInfo.makeData();
        mostRecent = (DataArithmetic)incomingDataInfo.makeData();
        mostRecentBlock = (DataArithmetic)incomingDataInfo.makeData();
        blockCorrelation = (DataArithmetic)incomingDataInfo.makeData();
        firstBlock = (DataArithmetic)incomingDataInfo.makeData();
        correlationSum = (DataArithmetic)incomingDataInfo.makeData();
        work = (DataArithmetic)incomingDataInfo.makeData();

        dataGroup = new DataGroup(new Data[] { mostRecent, average, error,
                        standardDeviation, mostRecentBlock, blockCorrelation});
        
        reset();
        
        DataInfoFactory factory = incomingDataInfo.getFactory();
        factory.getTags().add(getTag());
        String incomingLabel = incomingDataInfo.getLabel();
        factory.setLabel(incomingLabel+" most recent");
        DataInfo mostRecentInfo = factory.makeDataInfo();
        mostRecentInfo.addTag(mostRecentTag);
        factory.setLabel(incomingLabel+" avg");
        DataInfo averageInfo = factory.makeDataInfo();
        averageInfo.addTag(averageTag);
        factory.setLabel(incomingLabel+" error");
        DataInfo errorInfo = factory.makeDataInfo();
        errorInfo.addTag(errorTag);
        factory.setLabel(incomingLabel+" stddev");
        DataInfo standardDeviationInfo = factory.makeDataInfo();
        standardDeviationInfo.addTag(standardDeviationTag);
        factory.setLabel(incomingLabel+" most recent blk");
        DataInfo mostRecentBlockInfo = factory.makeDataInfo();
        mostRecentBlockInfo.addTag(mostRecentBlockTag);
        factory.setLabel(incomingLabel+" blk correlation");
        DataInfo correlationInfo = factory.makeDataInfo();
        correlationInfo.addTag(blockCorrelationTag);
        
        dataInfo = new DataInfoGroup(incomingLabel, incomingDataInfo.getDimension(), new DataInfo[]{
            mostRecentInfo, averageInfo, errorInfo, standardDeviationInfo, mostRecentBlockInfo,
            correlationInfo});
        dataInfo.addTag(getTag());
        return dataInfo;
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

    /**
     * Enumerated type that can be used to indicated the statistic to be taken
     * from the accumulator (e.g., average, error, current value, etc.). An
     * array of these types can be given to the addDataSink method to specify
     * the type of statistics to be given to the DataSink.
     */
    public static class StatType extends EnumeratedType {

        protected StatType(String label, int index) {
            super(label);
            this.index = index;
        }

        public static final StatType MOST_RECENT = new StatType("Latest value", 0);
        public static final StatType AVERAGE = new StatType("Average", 1);
        public static final StatType ERROR = new StatType("67% Confidence limits", 2);
        public static final StatType STANDARD_DEVIATION = new StatType("Standard deviation", 3);
        public static final StatType MOST_RECENT_BLOCK = new StatType("Latest block average", 4);
        public static final StatType BLOCK_CORRELATION = new StatType("Block correlation", 5);
        public static StatType[] choices() {
            return new StatType[] {MOST_RECENT,AVERAGE,ERROR,STANDARD_DEVIATION,MOST_RECENT_BLOCK,BLOCK_CORRELATION};
        }
        
        public final int index;
    }//end of ValueType


    protected DataArithmetic sum; //sum(blockSum/blkSize) = sum(blockAvg)
    protected DataArithmetic sumSquare;//sum(blockAvg^2)
    protected DataArithmetic blockSum;//block(value)
    protected DataArithmetic blockSumSq;//block(value^2)
    protected DataArithmetic sumSquareBlock;//sum(value^2)
    protected DataArithmetic mostRecent;//most recent value
    protected DataArithmetic mostRecentBlock;//most recent blockAvg
    protected DataArithmetic average, error, standardDeviation;
    protected DataArithmetic correlationSum, blockCorrelation, firstBlock;
    protected DataArithmetic work;
    protected DataGroup dataGroup;
    protected int count, blockCountDown;
    protected int blockSize;
    private final DataTag mostRecentTag, averageTag, errorTag, standardDeviationTag, mostRecentBlockTag, blockCorrelationTag;

}//end of AccumulatorAverage
