package etomica.util;

import java.io.Serializable;


/* History
 * 09/08/02 (DAK) added set/get methods for xMin, xMax, nBins
 * 08/04/04 (DAK,AJS,NRC) deleted DataSource.X methods; de-implemented DataSource.X.  Dimension-related material removed
 */

/**
 * Simple Histogram implementation with a static x range and number of bins.
 * If an x value is given that falls outside the histogram's x range, the 
 * value is dropped.  The x range and number of bins can be changed explicitly,
 * but doing so will reset the histogram (losing all previously collected data).
 */
public class HistogramSimple implements Histogram, java.io.Serializable {
	protected double deltaX;
	private int sum;
	protected int[] counts;
	protected double[] histogram;
    protected double xValues[];
    protected double xMin;
    protected double xMax;
    protected int nBins;

    /**
     * Makes a new histogram instance, with the range of x values given by
     * xRange.
     */
    public HistogramSimple(DoubleRange xRange) {
        this(100, xRange);
	}
	
    /**
     * Makes a new histogram instance, with the range of x values given by
     * xRange and n bins.
     */
    public HistogramSimple(int n, DoubleRange xRange) {
        nBins = n;
        counts = new int[n];
        histogram = new double[n];
        xValues = new double[n];
        setXRange(xRange);
    }
	
	public void reset() {
        //resets all histogram values and counts to zero
	    sum = 0;
	    for(int i=0; i<nBins; i++) {
	        counts[i] = 0;
	    }
	}
	
    public void addValue(double x) {
        //takes new value and updates histogram
        if(x >= xMin && x <= xMax) {
	        int i = (int)Math.floor(((x-xMin)/deltaX));
	        if(i == nBins){i--;}
            counts[i]++;
	    }
	    sum++;
    }
    
    public void setXRange(DoubleRange xRange) {
//        this.xRange = xRange;
        xMin = xRange.minimum();
        xMax = xRange.maximum();
        deltaX = (xMax-xMin)/nBins;
        reset();
    }
    
    public DoubleRange getXRange() {
        return new DoubleRange(xMin, xMax);
    }
      
    public int getNBins() {
        return nBins;
    }
    
    public void setNBins(int n) {
        this.nBins = n;
        counts = new int[n];
        histogram = new double[n];
        xValues = new double[n];
        deltaX = (xMax-xMin)/nBins;
        reset();
    }
    
    public double[] getHistogram() {
        //returns an array representing the present histogram
        if(sum != 0) {
		    for(int i=0; i<nBins; i++) {
		        histogram[i] = (double)counts[i]/((double)sum*deltaX);
		    }
        }
	    return histogram;
    }
    
    public int getCount() {
        return sum;
    }
 
    public double[] xValues() {
        return xValues;
    }
    
    public static class Factory implements Histogram.Factory, Serializable {
        public Factory(int n, DoubleRange xRange) {
            nBins = n;
            this.xRange = xRange;
        }
        
        public Histogram makeHistogram() {
            return new HistogramSimple(nBins, xRange);
        }
        
        private int nBins;
        private DoubleRange xRange;
    }
}
