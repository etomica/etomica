/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.types.DataDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.units.dimensions.Null;
import etomica.util.Statefull;

import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;

/**
 * Accumulator that keeps statistics for averaging and error analysis. The
 * statistics are incremented with every call to addData. Statistics kept are
 * <ul>
 * <li>the value last given to addData
 * <li>the mean, or average
 * <li>the statistical error (68% confidence limits), obtained by standard
 * error analysis of block averages
 * <li>the standard deviation
 * <li>the most recent block average
 * <li>the correlation between adjacent blocks (covariance divided by variance)
 * </ul>
 * Each statistic is recorded by a Data instance of the same type as the
 * incoming data. The output Data is a DataGroup formed from these Data
 * instances (in the order given above).
 * <p>
 * A <i>block </i> is a subset of the input data, formed by grouping (for
 * example) 100 or 1000 successive contributions to addData. Averages for each
 * block are considered to be independent of other block averages, and the
 * confidence interval for the overall average is obtained as the standard error
 * of the mean of these block averages.
 */
public abstract class AccumulatorAverage extends DataAccumulator implements Statefull {

    public static final StatType MOST_RECENT = new StatType("Latest value", 0);
    public static final StatType AVERAGE = new StatType("Average", 1);
    public static final StatType ERROR = new StatType("67% Confidence limits", 2);
    public static final StatType STANDARD_DEVIATION = new StatType("Standard deviation", 3);
    public static final StatType BLOCK_CORRELATION = new StatType("Block correlation", 4);
    private final DataTag mostRecentTag, averageTag, errorTag, standardDeviationTag, blockCorrelationTag;
    protected IData mostRecent;//most recent value
    protected IData average, error, standardDeviation;
    protected IData blockCorrelation;
    protected DataGroup dataGroup;
    protected long count, blockCountDown;
    protected long blockSize;
    protected boolean doIncludeACInError;

    /**
     * Default constructor sets block size to 1000 and sets the
     * interval for pushing the output data (pushInterval) to 100.
     */
    public AccumulatorAverage() {
        this(1000);
    }

    /**
     * Constructs with interval for pushing the output data set to 100.
     *
     * @param blockSize size of the block.
     */
    public AccumulatorAverage(long blockSize) {
        super();
        setBlockSize(blockSize);
        setPushInterval(100);
        mostRecentTag = new DataTag();
        averageTag = new DataTag();
        errorTag = new DataTag();
        standardDeviationTag = new DataTag();
        blockCorrelationTag = new DataTag();
    }

    /**
     * @return array of StatTypes identifying the statistics given in the DataGroup returned by getData().
     */
    public static StatType[] statChoices() {
        return new StatType[]{MOST_RECENT, AVERAGE, ERROR, STANDARD_DEVIATION, BLOCK_CORRELATION};
    }

    /**
     * @param statType the statistic for which the DataTag is requested.
     * @return the DataTag that identifies the specified statistic.
     */
    public DataTag getTag(StatType statType) {
        if (statType == MOST_RECENT) {
            return mostRecentTag;
        }
        if (statType == AVERAGE) {
            return averageTag;
        }
        if (statType == ERROR) {
            return errorTag;
        }
        if (statType == STANDARD_DEVIATION) {
            return standardDeviationTag;
        }
        if (statType == BLOCK_CORRELATION) {
            return blockCorrelationTag;
        }
        return null;
    }

    /**
     * Returns a specific element from the DataGroup returned
     * by getData().  Calling getData and extracting the data can be more
     * efficient, while this method can be more convenient.
     *
     * @return an IData object corresponding to the requested StatType.
     */
    public IData getData(StatType stat) {
        // populate data
        getData();
        if (dataGroup == null) return null;
        return dataGroup.getData(stat.index);
    }

    /**
     * Resets all sums to zero. All statistics are cleared.
     */
    public void reset() {
        count = 0;
        blockCountDown = blockSize;
        if (error == null) {
            return;
        }
        error.E(Double.NaN);
        mostRecent.E(Double.NaN);
        average.E(Double.NaN);
        standardDeviation.E(Double.NaN);
        blockCorrelation.E(Double.NaN);
    }

    /**
     * @return the current value of the block size.
     */
    public long getBlockSize() {
        return blockSize;
    }

    /**
     * Sets the size of the block used to group data for error analysis. Resets
     * statistics accumulated so far.
     *
     * @param newBlockSize new block size.
     */
    public void setBlockSize(long newBlockSize) {
        if (newBlockSize == 0) throw new IllegalArgumentException("block size must be positive");
        blockSize = newBlockSize;
        reset();
    }

    /**
     * Prepares the accumulator for input data.  Discards any previous
     * contributions to statistics.
     *
     * @param incomingDataInfo the DataInfo instance for the data that will be given to
     *                         addData
     * @return dataInfo object for the output of this accumulator.
     */
    public IDataInfo processDataInfo(IDataInfo incomingDataInfo) {
        standardDeviation = incomingDataInfo.makeData();
        average = incomingDataInfo.makeData();
        error = incomingDataInfo.makeData();
        mostRecent = incomingDataInfo.makeData();
        blockCorrelation = incomingDataInfo.makeData();

        dataGroup = new DataGroup(new IData[]{mostRecent, average, error,
                standardDeviation, blockCorrelation});

        reset();

        IDataInfoFactory factory = incomingDataInfo.getFactory();
        String incomingLabel = incomingDataInfo.getLabel();
        factory.setLabel(incomingLabel + " most recent");
        IDataInfo mostRecentInfo = factory.makeDataInfo();
        mostRecentInfo.addTag(mostRecentTag);
        factory.setLabel(incomingLabel + " avg");
        IDataInfo averageInfo = factory.makeDataInfo();
        averageInfo.addTag(averageTag);
        factory.setLabel(incomingLabel + " error");
        IDataInfo errorInfo = factory.makeDataInfo();
        errorInfo.addTag(errorTag);
        factory.setLabel(incomingLabel + " stddev");
        IDataInfo standardDeviationInfo = factory.makeDataInfo();
        standardDeviationInfo.addTag(standardDeviationTag);
        factory.setLabel(incomingLabel + " blk correlation");
        factory.setDimension(Null.DIMENSION);
        IDataInfo correlationInfo = factory.makeDataInfo();
        correlationInfo.addTag(blockCorrelationTag);

        dataInfo = new DataInfoGroup(incomingLabel, incomingDataInfo.getDimension(), new IDataInfo[]{
                mostRecentInfo, averageInfo, errorInfo, standardDeviationInfo, correlationInfo});
        dataInfo.addTags(incomingDataInfo.getTags());
        dataInfo.addTag(tag);
        return dataInfo;
    }

    /**
     * Adds a new DataSink that will receive a specific subset of the statistics
     * generated by this accumulator. The DataSink will be given a DataGroup
     * that has only the specified statistics in it. The output of this
     * accumulator is passed through a DataGroupFilter that removes the unwanted
     * Data.
     *
     * @param newDataSink the new DataSink
     * @param types       array indicating the statistics to be included in the
     *                    DataGroup sent to the sink.
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

    /**
     * @return the number of blocks that have been accumulated since construction or
     * the last call to reset.
     */
    public long getBlockCount() {
        return count;
    }

    /**
     * @return the number of samples that have been accumulated since construction or
     * the last call to reset.
     */
    public long getSampleCount() {
        return blockSize * count + (blockSize - blockCountDown);
    }

    /**
     * @return Flag indicating if auto-correlation is accounted for while calculating the uncertainty.
     * @see #setIncludeACInError
     */
    public boolean doIncludeACInError() {
        return doIncludeACInError;
    }

    /**
     * Include correction to the uncertainty from block correlation function as
     * prescribed in
     * Kolafa, Jiri(1986) 'Autocorrelations and subseries averages in Monte
     * Carlo Simulations', Molecular Physics (59) 1035
     *
     * @param newDoIncludeACInError auto-correlation is accounted if true, otherwise it is not.
     */
    public void setIncludeACInError(boolean newDoIncludeACInError) {
        doIncludeACInError = newDoIncludeACInError;
    }

    /**
     * Enumerated type that can be used to indicate the statistic to be taken
     * from the accumulator (e.g., average, error, current value, etc.). An
     * array of these types can be given to the addDataSink method to specify
     * the type of statistics to be given to the DataSink.
     */
    public static class StatType {

        public final int index;
        private final String label;

        protected StatType(String label, int index) {
            this.label = label;
            this.index = index;
        }

        @Override
        public String toString() {
            return this.label;
        }
    }

    public void saveState(Writer fw) throws IOException {
        fw.write(""+count+" "+blockCountDown+"\n");
        fw.write(""+mostRecent.getValue(0));
        for (int i=1; i<mostRecent.getLength(); i++) {
            fw.write(" "+mostRecent.getValue(i));
        }
        fw.write("\n");
    }

    public void restoreState(BufferedReader br) throws IOException {
        String[] bits = br.readLine().split(" ");
        count = Long.parseLong(bits[0]);
        blockCountDown = Long.parseLong(bits[1]);
        bits = br.readLine().split(" ");
        for (int i=0; i<mostRecent.getLength(); i++) {
            double x = Double.parseDouble(bits[0]);
            if (mostRecent instanceof DataDouble && i==0) {
                ((DataDouble) mostRecent).x = x;
            }
            else if (mostRecent instanceof DataDoubleArray) {
                ((DataDoubleArray)mostRecent).getData()[i] = x;
            }
            else {
                throw new RuntimeException("we only handle DataDouble and DataDoubleArray");
            }
        }
    }
}
