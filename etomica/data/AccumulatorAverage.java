/*
 * History
 * Created on Jul 26, 2004 by kofke
 */
package etomica.data;

import etomica.Constants;
import etomica.DataSink;
import etomica.DataTranslator;
import etomica.DataType;
import etomica.Default;
import etomica.units.Dimension;

/**
 * Accumulator that keeps statistics for averaging and error analysis.
 */
public class AccumulatorAverage extends DataAccumulator {

	public AccumulatorAverage() {
		super();
        allData = new double[numStats()][]; 
		setNData(0);
		setBlockSize(Default.BLOCK_SIZE);
        setDimension(Dimension.UNDEFINED);
        reset();
	}
	
	public void setBlockSize(int blockSize) {
		this.blockSize = blockSize;
        blockCountDown = blockSize;
	}
	public int getBlockSize() {
		return blockSize;
	}

    /**
     * Add the given values to the sums and block sums.
     */
    public void addData(double[] value) {
    	if(value.length != nData) setNData(value.length);
        for(int i=0; i<nData; i++) if(Double.isNaN(value[i])) return;
    	for(int i=nDataMinus1; i>=0; i--) {
            double v = value[i];
    		mostRecent[i] = v;
			blockSum[i] += v;
			blockSumSq[i] += v*v;
    	}
		if(--blockCountDown == 0) {//count down to zero to determine completion of block
		    doBlockSum();
        }
    }
    
    private void doBlockSum() {
        count++;
        blockCountDown = blockSize;
        for(int i=nDataMinus1; i>=0; i--) {             
            blockSum[i] /= blockSize;//compute block average
            sum[i] += blockSum[i];
            sumSquare[i] += blockSum[i]*blockSum[i];
            sumSquareBlock[i] += blockSumSq[i];
            //reset blocks
            mostRecentBlock[i] = blockSum[i];
            blockSum[i] = 0.0;
            blockSumSq[i] = 0.0;
        }
    }
    
    public double[] getData() {
        int currentBlockCount = blockSize - blockCountDown;
        double countFraction = (double)currentBlockCount/(double)blockSize;
        double currentCount = count + countFraction;
        if(count+currentBlockCount == 0) {
            setNaN(data);
        } else {
            for(int i=0; i<nData; i++) {
                double currentBlockAverage = blockSum[i]/currentBlockCount;
                double avg = (sum[i] + countFraction*currentBlockAverage)/currentCount;
                double avgSquared = avg*avg;
                double currentSumSquare = (sumSquare[i] + countFraction*currentBlockAverage*currentBlockAverage);
                double err = Math.sqrt(currentSumSquare/currentCount - avgSquared)/(currentCount-1);
                double stdev = Math.sqrt((sumSquareBlock[i]+blockSumSq[i])/(currentCount*blockSize) - avgSquared);
                double mrBlock = (!Double.isNaN(mostRecentBlock[i])) ? mostRecentBlock[i] : currentBlockAverage;
                data[0*nDataMinus1+i] = mostRecent[i];
                data[1*nDataMinus1+i] = average[i] = avg;
                data[2*nDataMinus1+i] = error[i] = err;
                data[3*nDataMinus1+i] = standardDeviation[i] = stdev;
                data[4*nDataMinus1+i] = mrBlock;
            }
        }
        return data;
    }
   
    protected void setNaN(double[] x) {
		for(int i=x.length-1; i>=0; i--) x[i] = Double.NaN;
	}
	protected void setZero(double[] x) {
		for(int i=x.length-1; i>=0; i--) x[i] = 0.0;
	}
	        
	/**
	 * Resets all sums to zero
	 */
    public void reset() {
        count = 0;
        setZero(sum);
        setZero(sumSquare);
        setZero(sumSquareBlock);
        setNaN(error);
        blockCountDown = blockSize;
        setZero(blockSum);
        setZero(blockSumSq);
        setNaN(mostRecent);
        setNaN(mostRecentBlock);
    }
    
    protected void setNData(int nData) {
    	this.nData = nData;
        if (nData > 1) {
            translator = new DataTranslatorArray(nData,numStats());
        }
        else {
            translator = DataTranslator.IDENTITY;
        }
        nDataMinus1 = nData-1;
        data = new double[numStats()*nData];
    	sum = redimension(nData, sum);
    	sumSquare = redimension(nData, sumSquare);
        sumSquareBlock = redimension(nData, sumSquareBlock);
    	standardDeviation = redimension(nData, standardDeviation);
    	average = redimension(nData, average);
    	error = redimension(nData, error);
    	blockSum = redimension(nData, blockSum);
        blockSumSq = redimension(nData, blockSumSq);
    	mostRecent = redimension(nData, mostRecent);
    	mostRecentBlock = redimension(nData, mostRecentBlock);
    	if(!saveOnRedimension) reset();
        allData[0] = mostRecent;
        allData[1] = average;
        allData[2] = error;
        allData[3] = standardDeviation;
        allData[4] = mostRecentBlock;
        for(int i=0; i<dataSinkList.length; i++) {
            if(dataSinkList[i] instanceof AccumulatorPusher) {
                ((AccumulatorPusher)dataSinkList[i]).setNData(nData);
            }
        }
    }
    
    /**
     * Creates a new array that is a redimensioning of the
     * given array, resized to the given integer value.
     * Truncates or pads with zeros as needed, and returns the
     * resized array.  Used by setNData.
     */
    protected double[] redimension(int n, double[] old) {
    	double[] newArray = new double[n];
    	if(saveOnRedimension && old != null) {
    		int k = (n > old.length) ? old.length : n;
            System.arraycopy(old, 0, newArray, 0, k);
    		//need to handle updating of counters, which should be different for new and old sums if saving on redimension
    		throw new etomica.exception.MethodNotImplementedException("Capability to save data on redimension not yet implemented"); 
    	}
    	return newArray;
    }

    public DataTranslator getTranslator() {
        return translator;
    }
    
	public DataType[] dataChoices() {return CHOICES;}
    
    public DataPusher makeDataPusher(Type[] types) {
       AccumulatorPusher newPusher = new AccumulatorPusher(types);
       addDataSink(newPusher.makeDataSink());
       return newPusher;
    }
    
    /**
	 * Typed constant that can be used to indicated the quantity
	 * to be taken from a meter (e.g., average, error, current value, etc.).
	 * Used primarily by Display objects.
	 */
	public static class Type extends etomica.DataType {
        protected Type(String label, int index) {
            super(label);
            this.index = index;
        }       
        public Constants.TypedConstant[] choices() {return CHOICES;}
        final int index;
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
    public static final Type MOST_RECENT_BLOCK = CHOICES[3];
    public static final Type STANDARD_DEVIATION = CHOICES[4];
	
    public int getCount() {
        return count;
    }

	/**
	 * @return Returns the saveOnRedimension.
	 */
	public boolean isSaveOnRedimension() {
		return saveOnRedimension;
	}
	/**
	 * @param saveOnRedimension The saveOnRedimension to set.
	 */
	public void setSaveOnRedimension(boolean saveOnRedimension) {
		this.saveOnRedimension = saveOnRedimension;
		if(saveOnRedimension) throw new IllegalArgumentException("Save on redimension not yet implemented correctly");
	}
    
    public int numStats() {
        return 5;
    }
	
    protected double[] sum, sumSquare, blockSum, blockSumSq, sumSquareBlock;
    protected double[] mostRecent;
    protected double[] mostRecentBlock;
    protected double[] average, error, standardDeviation;
    protected int count, blockCountDown;
    protected int blockSize;
    protected int nDataMinus1;
    protected int nData;
    protected boolean saveOnRedimension = false;
    
    //array concatenating mostRecent, average, etc. for return by getData
    protected double[] data;
    
    //the elements of allData point to mostRecent, average, etc. arrays (see setNData)
    protected final double[][] allData;
    
    protected DataTranslator translator;

    private class AccumulatorPusher extends DataPusher {
        
        AccumulatorPusher(Type[] types) {
            indexes = new int[types.length];
            selectedAllData = new double[types.length][];
            for(int i=0; i<indexes.length; i++) {
                indexes[i] = types[i].index;
            }
            setNData(AccumulatorAverage.this.nData);
        }

        //ignore argument, use data directly from outer class 
        protected void pushData() {
            if(nData == 1) {
                for(int i=0; i<indexes.length; i++) {
                    selectedData[i] = selectedAllData[i][0];
                }
            } else {
                for(int i=0, k=0; i<indexes.length; i++, k+=nData) {
                    System.arraycopy(selectedAllData[i], 0, selectedData, k, nData);
                }
            }
            pushData(selectedData);
        }
        
        void setNData(int nData) {
            selectedData = new double[indexes.length*nData];
            for(int i=0; i<indexes.length; i++) {
                selectedAllData[i] = allData[indexes[i]];
            }
            if (nData > 1) {
                selectedTranslator = new DataTranslatorArray(nData,indexes.length);
            }
            else {
                selectedTranslator = DataTranslator.IDENTITY;
            }
        }
        
        public DataTranslator getTranslator() {
            return selectedTranslator;
        }
        
        protected DataSink makeDataSink() {
            return new DataSink() {
                public void putData(double[] dummy) {
                    pushData();
                }
                public void setLabel(String s) {}
                public void setDimension(Dimension d) {}
            };
        }
        
        private double[] selectedData;
        private final double[][] selectedAllData;
        private final int[] indexes;
        private DataTranslator selectedTranslator;
    }//end of AccumulatorPipe
    
}//end of AccumulatorAverage
