/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.math.function.Function;
import etomica.math.function.IFunction;

import java.io.FileWriter;
import java.io.IOException;

/**
 * AccumulatorAverage that maintains a fixed block size.
 * <p>
 * This accumulator accepts any type of Data
 */
public class AccumulatorAverageFixedOutputFile extends AccumulatorAverage {

    /**
     * Default constructor sets block size to Default value, and sets the
     * interval for pushing the output data (pushInterval) to 100.
     */
    public AccumulatorAverageFixedOutputFile() {
        this(1000);
    }
    
    public AccumulatorAverageFixedOutputFile(int blockSize) {
        super(blockSize);

        // this chops negative values to be 0, used in cases where value
        // is non-zero due to roundoff and happens to be negative
        negativeChop = new IFunction(){
            public double f(double x) {
                return x > 0 ? x : 0;
            }
        };
    }

    /**
     * Add the given values to the sums and block sums. If any of the given data
     * values is NaN, method returns with no effect on accumulation sums.
     */
    public boolean addData(IData data) {
        if (data.isNaN())
            return false;

        mostRecent.E(data);
        currentBlockSum.PE(data);
        work.E(data);
        work.TE(data);
        sumSquare.PE(work);
        if (--blockCountDown == 0) {//count down to zero to determine
                                    // completion of block
            doBlockSum();
        }
        return true;
    }
    
    /**
     * Performs the block sum after <tt>blockSize</tt> calls to addData.
     */
    protected void doBlockSum() {
        blockCounter += blockSize;
    	count++;
        sum.PE(currentBlockSum);
        blockCountDown = blockSize;
        currentBlockSum.TE(1 / (double) blockSize);//compute block average
       
        if (filename!=null){
        	try{
        		fileWriter.write(blockCounter + " " +currentBlockSum.getValue(0)+"\n");
        	} catch (IOException e){
        	
        	}
        }
        work.E(currentBlockSum);
        work.TE(currentBlockSum);
        sumBlockSquare.PE(work);
        if (!mostRecentBlock.isNaN()) {
            work.E(currentBlockSum);
            work.TE(mostRecentBlock);
            correlationSum.PE(work);
        }
        else {
            firstBlock.E(currentBlockSum);
        }
        //reset blocks
        mostRecentBlock.E(currentBlockSum);
        currentBlockSum.E(0.0);
    }

    /**
     * Returns a DataGroup with Data instances holding the statistics kept by
     * this accumulator (as described in general comments for this class).
     */
    public IData getData() {
        if (sum == null)
            return null;
        if (count > 1) {
            // calculate block properties (these require 2 or more blocks)
            
            // block average (later discarded)
            average.E(sum);
            average.TE(1.0 / (count*blockSize));
            work.E(average);
            work.TE(average);
            error.E(sumBlockSquare);
            error.TE(1 / (double) count);
            error.ME(work);
            error.map(negativeChop);

            // error's intermediate value is useful for calculating block correlation
            blockCorrelation.E(average);
            blockCorrelation.TE(-2*count);
            blockCorrelation.PE(firstBlock);
            blockCorrelation.PE(mostRecentBlock);
            blockCorrelation.TE(average);
            blockCorrelation.PE(correlationSum);
            blockCorrelation.TE(1.0/(count-1));
            blockCorrelation.PE(work);
            blockCorrelation.DE(error);
            
            // ok, now finish up with error
            error.TE(1 / (double) (count - 1));
            error.map(Function.Sqrt.INSTANCE);
        }
        else {
            error.E(Double.NaN);
            blockCorrelation.E(Double.NaN);
        }
        
        long nTotalData = count*blockSize + (blockSize-blockCountDown);
        if (nTotalData > 0) {
            // now use *all* of the data
            average.E(sum);
            average.PE(currentBlockSum);
            average.TE(1.0 / nTotalData);
            work.E(average);
            work.TE(average);
            standardDeviation.E(sumSquare);
            standardDeviation.TE(1.0 / nTotalData);
            standardDeviation.ME(work);
            standardDeviation.map(negativeChop);
            standardDeviation.map(Function.Sqrt.INSTANCE);
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
        super.reset();
        if (sum == null) {
            return;
        }
        sum.E(0);
        sumBlockSquare.E(0);
        currentBlockSum.E(0);
        sumSquare.E(0);
        correlationSum.E(0);
        firstBlock.E(Double.NaN);
        mostRecentBlock.E(Double.NaN);
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
        sumBlockSquare = incomingDataInfo.makeData();
        currentBlockSum = incomingDataInfo.makeData();
        sumSquare = incomingDataInfo.makeData();
        firstBlock = incomingDataInfo.makeData();
        correlationSum = incomingDataInfo.makeData();
        mostRecentBlock = incomingDataInfo.makeData();
        work = incomingDataInfo.makeData();
        return super.processDataInfo(incomingDataInfo);
    }
    
    public void setFile(String fname){
    	this.filename = fname;
    	
    	try {
    		fileWriter = new FileWriter(fname);
    	} catch (IOException e){
    		fileWriter = null;
    	}
    }
    
    public void closeFile(){
    	try {
    		fileWriter.close();
    	} catch (IOException e){
    		
    	}
    }

    private static final long serialVersionUID = 1L;
    protected IData sum; //sum(value)
    protected IData sumBlockSquare;//sum(blockAvg^2)
    protected IData currentBlockSum;//block_sum(value)
    protected IData sumSquare;//sum(value^2)
    protected IData mostRecentBlock, correlationSum, firstBlock;
    protected IData work;
    protected final IFunction negativeChop;
    protected FileWriter fileWriter;
    private String filename;
    protected boolean toWriteFile = false;
    protected long blockCounter =0;

}
