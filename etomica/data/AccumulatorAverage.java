/*
 * History
 * Created on Jul 26, 2004 by kofke
 */
package etomica.data;

import etomica.Constants;
import etomica.Data;
import etomica.DataInfo;
import etomica.DataSink;
import etomica.Default;
import etomica.Constants.TypedConstant;
import etomica.data.types.DataArithmetic;
import etomica.data.types.DataGroup;
import etomica.utility.Function;

/**
 * Accumulator that keeps statistics for averaging and error analysis.
 */
public class AccumulatorAverage extends DataAccumulator {

	public AccumulatorAverage() {
		super();
		setBlockSize(Default.BLOCK_SIZE);
        setPushInterval(100);
	}
	
	public void setBlockSize(int blockSize) {
		this.blockSize = blockSize;
        blockCountDown = blockSize;
	}
	public int getBlockSize() {
		return blockSize;
	}
    
    public DataInfo getDataInfo() {
        return dataGroup.getDataInfo();
    }

    /**
     * Add the given values to the sums and block sums.  If any 
     * of the given data values is NaN, method returns with no 
     * effect on accumulation sums.
     */
    public void addData(Data data) {
        DataArithmetic value = (DataArithmetic)data;
        if(value.isNaN()) return;
        if (mostRecent == null) {
            initialize(data);
        }
  		mostRecent.E(data);
  	    blockSum.PE(value);
        work.E(data);
        work.TE(value);
        blockSumSq.PE(work);
		if(--blockCountDown == 0) {//count down to zero to determine completion of block
		    doBlockSum();
        }
    }
    
    protected void doBlockSum() {
        count++;
        blockCountDown = blockSize;
        blockSum.TE(1/(double)blockSize);//compute block average
        sum.PE(blockSum);
        work.E((Data)blockSum);
        work.TE(blockSum);
        sumSquare.PE(work);
        sumSquareBlock.PE(blockSumSq);
        //reset blocks
        mostRecentBlock.E((Data)blockSum);
        blockSum.E(0.0);
        blockSumSq.E(0.0);
    }
    
    public Data getData() {
        if (mostRecent == null) return null;
        int currentBlockCount = blockSize - blockCountDown;
        double countFraction = (double)currentBlockCount/(double)blockSize;
        double currentCount = count + countFraction;
        if(count+currentBlockCount > 0) {
//            double currentBlockAverage = blockSum[i]/currentBlockCount;
//            if (countFraction > 0) {
//                average = (sum[i] + countFraction*currentBlockAverage)/currentCount;
//            }
//            else {
            average.E((Data)sum);
            average.TE(1/(double)count);
            work.E((Data)average);
            work.TE(average);
            error.E((Data)sumSquare);
            error.TE(1/(double)count);
            error.ME(work);
            error.TE(1/(double)(count-1));
            error.map(Function.Sqrt.INSTANCE);
            standardDeviation.E((Data)sumSquareBlock);
            standardDeviation.PE(blockSumSq);
            standardDeviation.TE(1/currentCount*blockSize);
            standardDeviation.ME(work);
            standardDeviation.map(Function.Sqrt.INSTANCE);
//            mrBlock = (!Double.isNaN(mostRecentBlock[i])) ? mostRecentBlock[i] : currentBlockAverage;
        }
        return dataGroup;
    }
   
	/**
	 * Resets all sums to zero
	 */
    public void reset() {
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
    }
    
    protected void initialize(Data value) {
        sum = (DataArithmetic)value.makeCopy();
        sumSquare = (DataArithmetic)value.makeCopy();
        sumSquareBlock = (DataArithmetic)value.makeCopy();
        standardDeviation = (DataArithmetic)value.makeCopy();
        average = (DataArithmetic)value.makeCopy();
        error = (DataArithmetic)value.makeCopy();
        blockSum = (DataArithmetic)value.makeCopy();
        blockSumSq = (DataArithmetic)value.makeCopy();
        mostRecent = (DataArithmetic)value.makeCopy();
        mostRecentBlock = (DataArithmetic)value.makeCopy();
        work = (DataArithmetic)value.makeCopy();
        reset();
        dataGroup = new DataGroup(value.getDataInfo(),new Data[]{(Data)mostRecent,
                (Data)average,(Data)error,(Data)standardDeviation,(Data)mostRecentBlock});
    }
    
    public void addDataSink(DataSink dataSink, Type[] types) {
        int[] indexes = new int[types.length];
        for (int i=0; i<types.length; i++) {
            indexes[i] = types[i].index;
        }
        DataGroupFilter filter = new DataGroupFilter(indexes);
        addDataSink(filter);
        filter.addDataSink(dataSink);
    }
    
    public int getCount() {
        return count;
    }

    public TypedConstant[] dataChoices() {return CHOICES;}
    
    /**
	 * Typed constant that can be used to indicated the quantity
	 * to be taken from a meter (e.g., average, error, current value, etc.).
	 * Used primarily by Display objects.
	 */
	public static class Type extends TypedConstant {
        protected Type(String label, int index) {
            super(label);
            this.index = index;
        }       
        public Constants.TypedConstant[] choices() {return CHOICES;}
        public final int index;
    }//end of ValueType
    protected static final Type[] CHOICES = 
        new Type[] {
            new Type("Latest value", 0),
            new Type("Average", 1), 
            new Type("67% Confidence limits", 2),
            new Type("Standard deviation", 3),
            new Type("Latest block average", 4)};
    public static final Type MOST_RECENT = CHOICES[0];
    public static final Type AVERAGE = CHOICES[1];
    public static final Type ERROR = CHOICES[2];
    public static final Type STANDARD_DEVIATION = CHOICES[3];
    public static final Type MOST_RECENT_BLOCK = CHOICES[4];
	
    protected DataArithmetic sum, sumSquare, blockSum, blockSumSq, sumSquareBlock;
    protected DataArithmetic mostRecent;
    protected DataArithmetic mostRecentBlock;
    protected DataArithmetic average, error, standardDeviation;
    protected DataArithmetic work;
    protected DataGroup dataGroup;
    protected int count, blockCountDown;
    protected int blockSize;
    
}//end of AccumulatorAverage
