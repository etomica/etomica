/*
 * History
 * Created on Jul 26, 2004 by kofke
 */
package etomica.data;

import etomica.data.types.DataArithmetic;
import etomica.data.types.DataGroup;
import etomica.simulation.Simulation;
import etomica.units.Dimension;
import etomica.units.Undefined;
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
    }

    /**
     * Returns the DataInfo of the DataGroup that wraps the Data statistics.
     */
    public DataInfo getDataInfo() {
        return dataGroup.getDataInfo();
    }

    /**
     * Checks that incoming Data implements DataArithmetic, and returns null if
     * this is so. Otherwise throws a ClassCastException, as there is no data
     * caster to DataArithmetic.
     */
    public DataProcessor getDataCaster(DataInfo incomingDataInfo) {
        if (DataArithmetic.class.isAssignableFrom(incomingDataInfo.getDataClass())) {
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
        work.E((Data) blockSum);
        work.TE(blockSum);
        sumSquare.PE(work);
        if (!mostRecentBlock.isNaN()) {
            work.E((Data)blockSum);
            work.TE(mostRecentBlock);
            correlationSum.PE(work);
        }
        else {
            firstBlock.E((Data)blockSum);
        }
        sumSquareBlock.PE(blockSumSq);
        //reset blocks
        mostRecentBlock.E((Data) blockSum);
        blockSum.E(0.0);
        blockSumSq.E(0.0);
    }

    /**
     * Returns a DataGroup with Data instances holding the statistics kept by
     * this accumulator (as described in general comments for this class).
     */
    public Data getData() {
        if (mostRecent == null)
            return null;
        //        int currentBlockCount = blockSize - blockCountDown;
        //        double countFraction = (double)currentBlockCount/(double)blockSize;
        //        double currentCount = count + countFraction;
        if (count > 0) {
            //            double currentBlockAverage = blockSum[i]/currentBlockCount;
            //            if (countFraction > 0) {
            //                average = (sum[i] +
            // countFraction*currentBlockAverage)/currentCount;
            //            }
            //            else {
            average.E((Data) sum);
            average.TE(1 / (double) count);
            work.E((Data) average);
            work.TE(average);
            error.E((Data) sumSquare);
            error.TE(1 / (double) count);
            error.ME(work);

            // error's intermediate value is useful for calculating block correlation
            blockCorrelation.E((Data)sum);
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
            standardDeviation.E((Data) sumSquareBlock);
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
        DataFactory factory = incomingDataInfo.getDataFactory();

        Dimension dimSquared = Undefined.DIMENSION;//can change this when units
                                                   // facility isbetter
                                                   // developed
        sum = (DataArithmetic) factory.makeData("blkAvg sum", incomingDataInfo.getDimension());
        sumSquare = (DataArithmetic) factory.makeData("blkAvgSqr sum", dimSquared);
        sumSquareBlock = (DataArithmetic) factory.makeData("sum value^2", incomingDataInfo.getDimension());
        standardDeviation = (DataArithmetic) factory.makeData("stddev", incomingDataInfo.getDimension());
        average = (DataArithmetic) factory.makeData("avg", incomingDataInfo.getDimension());
        error = (DataArithmetic) factory.makeData("error", incomingDataInfo.getDimension());
        blockSum = (DataArithmetic) factory.makeData("blk value", incomingDataInfo.getDimension());
        blockSumSq = (DataArithmetic) factory.makeData("blk value^2", dimSquared);
        mostRecent = (DataArithmetic) factory.makeData("most recent", incomingDataInfo.getDimension());
        mostRecentBlock = (DataArithmetic) factory.makeData("most recent blk", incomingDataInfo.getDimension());
        blockCorrelation = (DataArithmetic) factory.makeData("blk correlation", incomingDataInfo.getDimension());
        firstBlock = (DataArithmetic) factory.makeData("first blk", incomingDataInfo.getDimension());
        correlationSum = (DataArithmetic) factory.makeData("correlation sum", incomingDataInfo.getDimension());
        work = (DataArithmetic) factory.makeData("scratch", Undefined.DIMENSION);

        reset();
        dataGroup = new DataGroup(incomingDataInfo.getLabel(),
                new Data[] { (Data) mostRecent, (Data) average, (Data) error,
                        (Data) standardDeviation, (Data) mostRecentBlock,(Data)blockCorrelation});
        return dataGroup.getDataInfo();
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

}//end of AccumulatorAverage
