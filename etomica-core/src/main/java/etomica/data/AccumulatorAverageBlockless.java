/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.math.function.IFunction;
import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.util.EnumeratedType;
import etomica.math.function.Function;

/**
 * AccumulatorAverage that does not use blocks (best used when each sample
 * is independent).  This accumulator tracks most recent, average and standard
 * deviation.
 * <p>
 * This accumulator accepts any type of Data
 */
public class AccumulatorAverageBlockless extends DataAccumulator {

    /**
     * Default constructor sets block size to Default value, and sets the
     * interval for pushing the output data (pushInterval) to 100.
     */
    public AccumulatorAverageBlockless() {

        setPushInterval(100);
        mostRecentTag = new DataTag();
        averageTag = new DataTag();
        standardDeviationTag = new DataTag();

        // this chops negative values to be 0, used in cases where value
        // is non-zero due to roundoff and happens to be negative
        negativeChop = new IFunction(){
            public double f(double x) {
                return !(x <= 0) ? x : 0;
            }
        };
    }
    
    public DataTag getTag(StatType statType) {
        if (statType == MOST_RECENT) {
            return mostRecentTag;
        }
        if (statType == AVERAGE) {
            return averageTag;
        }
        if (statType == STANDARD_DEVIATION) {
            return standardDeviationTag;
        }
        return null;
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
    public boolean addData(IData data) {
        if (data.isNaN())
            return false;

        mostRecent.E(data);
        sum.PE(data);
        work.E(data);
        work.TE(data);
        sumSquare.PE(work);
        count++;
        return true;
    }

    /**
     * Returns a DataGroup with Data instances holding the statistics kept by
     * this accumulator (as described in general comments for this class).
     */
    public IData getData() {
        if (sum == null)
            return null;

        if (count > 0) {
            average.E(sum);
            average.DE(count);
            work.E(average);
            work.TE(average);
            if (count > 1) {
                standardDeviation.E(sumSquare);
                standardDeviation.DE(count);
                standardDeviation.ME(work);
                standardDeviation.map(negativeChop);
                standardDeviation.map(Function.Sqrt.INSTANCE);
            }
            else {
                standardDeviation.E(Double.NaN);
            }
        }
        else {
            average.E(Double.NaN);
            standardDeviation.E(Double.NaN);
        }
        return dataGroup;
    }

    /**
     * Resets all sums to zero. All statistics are cleared.
     */
    public void reset() {
        /**
         * Resets all sums to zero. All statistics are cleared.
         */
        count = 0;
        if (sum== null) {
            return;
        }
        mostRecent.E(Double.NaN);
        average.E(Double.NaN);
        standardDeviation.E(Double.NaN);

        sum.E(0);
        sumSquare.E(0);
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
        work = incomingDataInfo.makeData();
        standardDeviation = incomingDataInfo.makeData();
        average = incomingDataInfo.makeData();
        mostRecent = incomingDataInfo.makeData();

        dataGroup = new DataGroup(new IData[] { mostRecent, average, standardDeviation});
        
        reset();
        
        IDataInfoFactory factory = incomingDataInfo.getFactory();
        String incomingLabel = incomingDataInfo.getLabel();
        factory.setLabel(incomingLabel+" most recent");
        IDataInfo mostRecentInfo = factory.makeDataInfo();
        mostRecentInfo.addTag(mostRecentTag);
        factory.setLabel(incomingLabel+" avg");
        IDataInfo averageInfo = factory.makeDataInfo();
        averageInfo.addTag(averageTag);
        factory.setLabel(incomingLabel+" stddev");
        IDataInfo standardDeviationInfo = factory.makeDataInfo();
        standardDeviationInfo.addTag(standardDeviationTag);
        
        dataInfo = new DataInfoGroup(incomingLabel, incomingDataInfo.getDimension(), new IDataInfo[]{
            mostRecentInfo, averageInfo, standardDeviationInfo});
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

    public long getSampleCount() {
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

        private static final long serialVersionUID = 1L;
        public final int index;
    }

    public static final StatType MOST_RECENT = new StatType("Latest value", 0);
    public static final StatType AVERAGE = new StatType("Average", 1);
    public static final StatType STANDARD_DEVIATION = new StatType("Standard deviation", 2);
    public static StatType[] statChoices() {
        return new StatType[] {MOST_RECENT,AVERAGE,STANDARD_DEVIATION};
    }

    public AllData getRawData() {
        return new AllData(count, sum, sumSquare);
    }

    public void setRawData(AllData ad) {
        count = ad.count;
        sum.E(ad.sum);
        sumSquare.E(ad.sumSquare);
    }

    public static class AllData {
        public final IData sum, sumSquare;
        public final long count;
        public AllData(long count, IData sum, IData sumSquare) {
            this.sum = sum;
            this.sumSquare = sumSquare;
            this.count = count;
        }
    }

    protected IData mostRecent;//most recent value
    protected IData average, standardDeviation;
    protected DataGroup dataGroup;
    protected long count;
    private final DataTag mostRecentTag, averageTag, standardDeviationTag;

    protected IData sum; //sum(value)
    protected IData sumSquare;//sum(value^2)
    protected IData work;
    protected final IFunction negativeChop;
}
