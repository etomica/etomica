/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.math.function.IFunction;
import etomica.data.AccumulatorAverageBlockless.AllData;
import etomica.data.types.DataGroup;
import etomica.util.EnumeratedType;
import etomica.math.function.Function;

/**
 * AccumulatorAverage that does not use blocks (best used when each sample
 * is independent).  This accumulator tracks average and standard deviation.
 * The returned data is created on the fly in order to save on memory.
 * <p>
 * This accumulator accepts any type of Data
 */
public class AccumulatorAverageBlocklessSlim extends DataAccumulator {
    
    static {
        // this chops negative values to be 0, used in cases where value
        // is non-zero due to roundoff and happens to be negative
        negativeChop = new IFunction(){
            public double f(double x) {
                return !(x <= 0) ? x : 0;
            }
        };
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

        IData average = inDataInfo.makeData();
        IData standardDeviation = inDataInfo.makeData();
        DataGroup dataGroup = new DataGroup(new IData[]{average,standardDeviation});
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

        sum.E(0);
        sumSquare.E(0);
    }

    /**
     * Prepares the accumulator for input data.  Discards any previous 
     * contributions to statistics.
     * 
     * @param inputDataInfo
     *            the DataInfo instance for the data that will be given to
     *            addData
     */
    public void putDataInfo(IDataInfo inputDataInfo) {
        inDataInfo = inputDataInfo;
        sum = inputDataInfo.makeData();
        sumSquare = inputDataInfo.makeData();
        work = inputDataInfo.makeData();
        
        reset();
    }

    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        return null;
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

    public static final StatType AVERAGE = new StatType("Average", 0);
    public static final StatType STANDARD_DEVIATION = new StatType("Standard deviation", 1);
    public static StatType[] statChoices() {
        return new StatType[] {AVERAGE,STANDARD_DEVIATION};
    }

    public AllData getRawData() {
        return new AllData(count, sum, sumSquare);
    }

    public void setRawData(AllData ad) {
        count = ad.count;
        sum.E(ad.sum);
        sumSquare.E(ad.sumSquare);
    }

    protected IDataInfo inDataInfo;
    protected long count;

    protected IData sum; //sum(value)
    protected IData sumSquare;//sum(value^2)
    protected IData work;// used to compute value^2
    protected static final IFunction negativeChop;
}
