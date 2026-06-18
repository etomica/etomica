/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.histogram;


import etomica.math.DoubleRange;

/**
 * Histogram implementation with dynamic x range and number of bins.  If an x 
 * value is given that falls outside the histogram's x range, the x range 
 * expands to accomodate any values outside the current range maintaining the 
 * bin width at a fixed value.  The number of bins and the X range can be set 
 * explicitly, but only if the change would not cause old data to be dropped.  
 * Unless one is given, there is no initial range.  The first incoming piece of
 * data sets the intial range.
 *
 * @author Andrew Schultz
 */
public class HistogramExpanding extends HistogramSimple {

    /**
     * Default constructor, making a histogram with deltaX = 1.0.
     */
    public HistogramExpanding() {this(1.0);}

    /**
     * Constructs a new instance with the give deltaX bin width and no range.
     */
    public HistogramExpanding(double deltaX) {
        this(deltaX, new DoubleRange(Double.NaN, Double.NaN));
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
	
    public void reset() {
        // undefine the range so that the next data point redefines the range.
        setXRange(new DoubleRange(Double.NaN, Double.NaN));
    }
    
    public void addValue(double x) {     //takes new value and updates histogram
    	if (Double.isNaN(xMin)) {
    		// the first piece of data, so create the initial range to match 
    		// this point
    		double newXMin = (int)(x/deltaX)*deltaX;
    		if (x < 0) { 
    			newXMin -= deltaX;
    		}
    		setXRange(new DoubleRange(newXMin,newXMin+deltaX));
    	}
    	else if (x < xMin) {
    	    if(x>xMax)throw new RuntimeException(x + " " + xMin + " " + xMax);
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
        if (Double.isNaN(newXMin)) {
        	// we're setting the range to nothing and then we'll wait for
        	// data to come in and use that to initialize the range
        	xMin = Double.NaN;
        	xMax = Double.NaN;
        	counts = new long[0];
            histogram = new double[0];
            xValues = new double[0];
        	nBins = 0;
        	return;
        }
        newXMin = Math.floor(newXMin/deltaX) * deltaX;
        newXMax = Math.ceil(newXMax/deltaX) * deltaX;
        int newNBins = (int)Math.round((newXMax - newXMin) / deltaX);
        if(newNBins <= 0) throw new RuntimeException(newNBins + " "+ newXMin + " "+ newXMax + " "+ xRange);
        // if called from the constructor, we have no data to copy, so skip this
        if (xMax != 0 || xMin != 0) {
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
            long[] newCounts = new long[newNBins];
            for (int i=Math.max(0,minDiff); i<nBins-Math.max(0,maxDiff); i++) {
                newCounts[i-minDiff] = counts[i];
            }
            counts = newCounts;
        }
        else {
            counts = new long[newNBins];
        }
        xMin = newXMin;
        xMax = newXMax;
        nBins = newNBins;
        histogram = new double[nBins];
        xValues = new double[nBins];

        for(int i=0; i<nBins; i++) {
            xValues[i] = xMin + deltaX * (i+0.5);
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

    public double getDeltaX(){
    	return deltaX;
    }
    
    public void setDeltaX(double dx){
    	if (dx == deltaX) return;
    	deltaX = dx;
    	
    	setXRange(new DoubleRange(xMin, xMax));
    }
    
    private static final long serialVersionUID = 1L;
}
