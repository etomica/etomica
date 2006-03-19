package etomica.util;

import java.io.Serializable;


/**
 * Histogram implementation with dynamic x range and number of bins.  If an x 
 * value is given that falls outside the histogram's x range, the x range 
 * expands to accomodate any values outside the current range maintaining the 
 * bin width at a fixed value.  The number of bins and the X range can be set 
 * explicitly, but only if the change would not cause old data to be dropped.  
 *
 * @author Andrew Schultz
 */
public class HistogramExpanding extends HistogramSimple {

    /**
     * Default constructor, making a histogram with deltaX = 1.0.
     */
    public HistogramExpanding() {this(1.0);}

    /**
     * Constructs a new instance with the give deltaX bin width.
     */
    public HistogramExpanding(double deltaX) {
        this(deltaX, new DoubleRange(0,10));
    }
    
    /**
     * Constructs a new instance with the given deltaX bin width and
     * range of X values.
     */
    public HistogramExpanding(double deltaX, DoubleRange xRange) {
        // super constructor here is useless, but we have to call it
        super(10, xRange);
        // we'd like deltaX to be final, but it's not.  We override both 
        // methods that might change it
        this.deltaX = deltaX;
        setXRange(xRange);
    }

    /**
     * Makes a new histogram instance, with n bins covering the
     * given range of X values.
     */
    public HistogramExpanding(int n, DoubleRange xRange) {
        this((xRange.maximum() - xRange.minimum()) / n, xRange);
    }
	
    public void addValue(double x) {     //takes new value and updates histogram
        if (x < xMin) {
            setXRange(new DoubleRange(x,xMax));
        }
        else if(x > xMax) {
            setXRange(new DoubleRange(xMin,x));
        }

        super.addValue(x);
    }
    
    public void setXRange(DoubleRange xRange) {
        if (deltaX == 0) {
            //called from super constructor
            return;
        }
        double newXMin = xRange.minimum();
        double newXMax = xRange.maximum();
        newXMin = Math.floor(newXMin/deltaX) * deltaX;
        newXMax = Math.ceil(newXMax/deltaX) * deltaX;
        int newNBins = (int)Math.round((newXMax - newXMin) / deltaX);
        // how many bins would get dropped on the left side
        int minDiff = (int)Math.round((newXMin - xMin) / deltaX);
        for (int i=0; i<minDiff; i++) {
            // loop over bins on the left getting dropped
            if (counts[i] != 0) {
                throw new IllegalStateException("complying would drop data");
            }
        }
        // how many bins would get dropped on the right side
        int maxDiff = (int)Math.round((xMax - newXMax) / deltaX);
        for (int i=nBins-1; i>nBins-maxDiff-1; i--) {
            // loop over bins on the right getting dropped
            if (counts[i] != 0) {
                throw new IllegalStateException("complying would drop data");
            }
        }
        //copy old data
        int[] newCounts = new int[newNBins];
        for (int i=Math.max(0,minDiff); i<nBins-Math.max(0,maxDiff); i++) {
            newCounts[i-minDiff] = counts[i];
        }
        counts = newCounts;
        xMin = newXMin;
        xMax = newXMax;
        nBins = newNBins;
        histogram = new double[nBins];
        xValues = new double[nBins];
        for(int i=0; i<nBins; i++) {
            xValues[i] = xMin + (i+0.5)*deltaX;
        }
    }

    /**
     * Sets the number of bins to n.  Adds or removes bins from the end of 
     * the X range.  For more flexibility, call setXRange directly.   
     */
    public void setNBins(int n) {
        if (n == nBins) return;
        // subtract a half to prevent roundup to n+1
        // setXRange will the do the "Right Thing" (round up to n)
        setXRange(new DoubleRange(xMin,xMin+deltaX*(n-0.5)));
    }

    public static final Histogram.Factory FACTORY = new Factory(1.0);
    
    public static class Factory implements Histogram.Factory, Serializable {
        public Factory(double deltaX) {
            this.deltaX = deltaX;
        }
        
        public Histogram makeHistogram() {
            return new HistogramExpanding(deltaX);
        }
        
        private final double deltaX;
    }

}
