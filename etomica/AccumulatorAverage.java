/*
 * History
 * Created on Jul 26, 2004 by kofke
 */
package etomica;

/**
 * @author kofke
 *
 * Accumulator that keeps statistics for averaging and error analysis.
 */
public class AccumulatorAverage extends Accumulator implements DataSourceMultitype {

	/**
	 * @param parentElement
	 * @param dataSource
	 */
	public AccumulatorAverage() {
		super();
		reset();
	}

    /**
     * Add the given values to the sums and block sums
     */
    public void add(double[] value) {
    	if(value.length != nData) setNData(value.length);
    	for(int i=nDataMinus1; i>=0; i--) {       		
    		mostRecent[i] = value[i];
//				if(Double.isNaN(value[i])) continue;				
			blockSum[i] += value[i];
			
    	}
		if(--blockCountDown == 0) {//count down to zero to determine completion of block
	       	for(int i=nDataMinus1; i>=0; i--) {       		
				blockSum[i] /= blockSize;//compute block average
			    sum[i] += blockSum[i];
			    sumSquare[i] += blockSum[i]*blockSum[i];
	            //reset blocks
	            mostRecentBlock[i] = blockSum[i];
	            blockSum[i] = 0.0;
			}
		    count++;
            blockCountDown = blockSize;
        }
    }
    
    /**
     * Compute and return the average from the present value of the accumulation sum
     */
    public double[] average() {
        int blockCount = blockSize - blockCountDown;
        if(count+blockCount > 0) setNaN(average);
        else {
	       	for(int i=nDataMinus1; i>=0; i--) {       		
	        	average[i] = sum[i]/(double)count;
//		            average[i] = (sum[i] + blockSum[i]/blockSize)/(double)(count + (double)blockCount/(double)blockSize) 
	        }
        }
        return average;
    }
	
	/**
	 * Return the 67% confidence limits of the average based on variance in block averages.
	 */
    public double[] error() {
    	if(count > 1) {
           	for(int i=nDataMinus1; i>=0; i--) {       		
		        double avg = sum[i]/(double)count;	    		
		        error[i] = Math.sqrt((sumSquare[i]/(double)count - avg*avg)/(double)(count-1));	    		
	    	}
    	} else {
    		setNaN(error);
    	}
        return error;
    }
    
    /**
     * Compute and return the standard deviation of the recorded data.
     */
     public double[] standardDeviation() {
     	if(count > 0) {
           	for(int i=nDataMinus1; i>=0; i--) {       		
           		double avg = sum[i]/(double)count;
           		standardDeviation[i] = Math.sqrt(sumSquare[i]/(double)count - avg*avg);
           	}
     	} else setNaN(standardDeviation);
     	return standardDeviation;
     }
    
    /**
     * Returns the value last passed to the add method
     */
    public double[] mostRecent() {return mostRecent;}
	
	/**
	 * Returns the value of the most recent block average.
	 */
	public double[] mostRecentBlock() {
       	for(int i=nDataMinus1; i>=0; i--) {       		
       		mostRecentBlock[i] = (mostRecentBlock[i]!=Double.NaN) ? mostRecentBlock[i] : blockSum[i]/(blockSize - blockCountDown);
       	}
       	return mostRecentBlock;
    }
	
	private void setNaN(double[] x) {
		for(int i=x.length-1; i>=0; i--) x[i] = Double.NaN;
	}
	private void setZero(double[] x) {
		for(int i=x.length-1; i>=0; i--) x[i] = 0.0;
	}
	        
	/**
	 * Resets all sums to zero
	 */
    public void reset() {
        count = 0;
        setZero(sum);
        setZero(sumSquare);
        setNaN(error);
        blockCountDown = blockSize;
        setZero(blockSum);
        setNaN(mostRecent);
        setNaN(mostRecentBlock);
    }
    
    private void setNData(int nData) {
    	this.nData = nData;
    	nDataMinus1 = nData-1;
    	sum = redimension(nData, sum);
    	sumSquare = redimension(nData, sumSquare);
    	average = redimension(nData, average);
    	error = redimension(nData, error);
    	mostRecent = redimension(nData, mostRecent);
    	mostRecentBlock = redimension(nData, mostRecentBlock);
    	if(!saveOnRedimension) reset();
    }
    
    /**
     * Creates a new array that is a redimensioning of the
     * given array, resized to the given integer value.
     * Truncates or pads with zeros as needed, and returns the
     * resized array.  Used by setNData.
     */
    private double[] redimension(int n, double[] old) {
    	double[] newArray = new double[n];
    	if(saveOnRedimension) {
    		int k = (n > old.length) ? old.length : n;
    		for(int i=0; i<k; i++) {
    			newArray[i] = old[i];
    		}
    		//need to handle updating of counters, which should be different for new and old sums if saving on redimension
    		throw new etomica.exception.MethodNotImplementedException("Capability to save data on redimension not yet implemented"); 
    	}
    	return newArray;
    }
        	 
    /**
    * Accessor method to indicate if the meter should keep a histogram of all measured values.
    * Default is false (do not keep histogram).
    * Does not take an argument; true/falue value is obtained from outer-class meter
    */
    //need way to update name if meter name changes
    public void setHistogramming() {
        if(histogramming && histogram == null) {
            histogram = new HistogramSimple();
            histogram.setName(MeterAbstract.this.toString() + ":HistogramSimple");
            histogram.setLabel(MeterAbstract.this.getLabel() + " histogram");
            histogram.setXLabel(MeterAbstract.this.getLabel());
            histogram.setXDimension(MeterAbstract.this.getDimension());
        }
    }
    
   /**
    * Accessor method to indicate if the meter should keep a history of all measured values.
    * Default is false (do not keep history).
    * Does not take an argument; true/falue value is obtained from outer-class meter
    */
    public void setHistorying() {
        if(historying && history == null) {
            history = new History(Default.HISTORY_PERIOD);
            history.setName(MeterAbstract.this.toString() + ":History");
            history.setLabel(MeterAbstract.this.getLabel() + " history");
        }
    }

	public double[] getData() {
	    return getData(AVERAGE);
	}
	
	public double[] getData(DataType type) {
		return getData((Type)type);
	}
	/**
	 * Returns the value indicated by the argument.
	 */
	public double[] getData(Type type) {
	    if(type==AVERAGE) return average();
	    if(type==MOST_RECENT) return mostRecent();
	    if(type==MOST_RECENT_BLOCK) return mostRecentBlock();
	    if(type==ERROR) return error();
	    if(type==STANDARD_DEVIATION) return standardDeviation();
		/*if(type==null)*/ return getData();
	}
	
	public DataType[] dataChoices() {return CHOICES;}
	
	/**
	 * Typed constant that can be used to indicated the quantity
	 * to be taken from a meter (e.g., average, error, current value, etc.).
	 * Used primarily by Display objects.
	 */
	public static class Type extends etomica.DataType {
        private Type(String label) {super(label);}       
        public Constants.TypedConstant[] choices() {return (Constants.TypedConstant[])CHOICES;}
    }//end of ValueType
    private static final Type[] CHOICES = 
        new Type[] {
            new Type("Average"), new Type("67% Confidence limits"),
            new Type("Latest value"),
            new Type("Latest block average"), new Type("Standard deviation")};
    public static final Type AVERAGE = CHOICES[0];
    public static final Type ERROR = CHOICES[1];
    public static final Type MOST_RECENT = CHOICES[2];
    public static final Type MOST_RECENT_BLOCK = CHOICES[3];
    public static final Type STANDARD_DEVIATION = CHOICES[4];
	

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
	
    private double[] sum, sumSquare, blockSum;
    private double[] mostRecent;
    private double[] mostRecentBlock;
    private double[] average, error, standardDeviation;
    private int count, blockCountDown;
    private int blockSize;
    private int nDataMinus1;
    private boolean saveOnRedimension = false;

}
