package etomica.data;

import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.util.EnumeratedType;

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
 */
public abstract class AccumulatorAverage extends DataAccumulator {

    /**
     * Default constructor sets block size to Default value, and sets the
     * interval for pushing the output data (pushInterval) to 100.
     */
    public AccumulatorAverage() {
        this(1000);
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

    public DataTag getTag(StatType statType) {
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
        if (statType == StatType.BLOCK_CORRELATION) {
            return blockCorrelationTag;
        }
        return null;
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
        blockCountDown = blockSize;
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
        standardDeviation = incomingDataInfo.makeData();
        average = incomingDataInfo.makeData();
        error = incomingDataInfo.makeData();
        mostRecent = incomingDataInfo.makeData();
        blockCorrelation = incomingDataInfo.makeData();

        dataGroup = new DataGroup(new Data[] { mostRecent, average, error,
                        standardDeviation, blockCorrelation});
        
        reset();
        
        IDataInfoFactory factory = incomingDataInfo.getFactory();
        String incomingLabel = incomingDataInfo.getLabel();
        factory.setLabel(incomingLabel+" most recent");
        IDataInfo mostRecentInfo = factory.makeDataInfo();
        mostRecentInfo.addTag(mostRecentTag);
        factory.setLabel(incomingLabel+" avg");
        IDataInfo averageInfo = factory.makeDataInfo();
        averageInfo.addTag(averageTag);
        factory.setLabel(incomingLabel+" error");
        IDataInfo errorInfo = factory.makeDataInfo();
        errorInfo.addTag(errorTag);
        factory.setLabel(incomingLabel+" stddev");
        IDataInfo standardDeviationInfo = factory.makeDataInfo();
        standardDeviationInfo.addTag(standardDeviationTag);
        factory.setLabel(incomingLabel+" most recent blk");
        IDataInfo mostRecentBlockInfo = factory.makeDataInfo();
        mostRecentBlockInfo.addTag(mostRecentBlockTag);
        factory.setLabel(incomingLabel+" blk correlation");
        IDataInfo correlationInfo = factory.makeDataInfo();
        correlationInfo.addTag(blockCorrelationTag);
        
        dataInfo = new DataInfoGroup(incomingLabel, incomingDataInfo.getDimension(), new IDataInfo[]{
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
        public static final StatType BLOCK_CORRELATION = new StatType("Block correlation", 4);
        public static StatType[] choices() {
            return new StatType[] {MOST_RECENT,AVERAGE,ERROR,STANDARD_DEVIATION,BLOCK_CORRELATION};
        }
        
        private static final long serialVersionUID = 1L;
        public final int index;
    }


    private static final long serialVersionUID = 1L;
    protected Data mostRecent;//most recent value
    protected Data average, error, standardDeviation;
    protected Data blockCorrelation;
    protected DataGroup dataGroup;
    protected int count, blockCountDown;
    protected int blockSize;
    private final DataTag mostRecentTag, averageTag, errorTag, standardDeviationTag, mostRecentBlockTag, blockCorrelationTag;

}
