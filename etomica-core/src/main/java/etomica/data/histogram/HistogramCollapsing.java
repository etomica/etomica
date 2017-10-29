/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.histogram;

import etomica.math.DoubleRange;

/**
 * Simple Histogram implementation with a static number of bins, but dynamic
 * x range. If an x value is given that falls outside the histogram's x range, 
 * the histogram's range is expanded to include the given x value.  The 
 * previously collected histogram is collapsed into fewer bins.  The x range 
 * and number of bins can be changed explicitly, but doing so will reset the 
 * histogram (losing all previously collected data).
 * 
 * @author Andrew Schultz
 */
public class HistogramCollapsing extends HistogramSimple {
    protected double firstValue;
    protected long firstValueCount;

    /**
     * Default constructor, making a 100-bin histogram.
     */
    public HistogramCollapsing() {
        this(100);
    }

    /**
     * Makes a new histogram instance having a number of bins given by the argument.
     */
    public HistogramCollapsing(int n) {
        // range is bogus, but we'll drop it anyway.
        super(n, new DoubleRange(0,0));
        firstValue = Double.NaN;
        firstValueCount = 0;
    }

	public void reset() {                //resets all histogram values and counts to zero
	    super.reset();
	    firstValue = Double.NaN;
	    firstValueCount = 0;
	}
	
    public void addValue(double x) {     //takes new value and updates histogram
        if (Double.isNaN(x)) return;
        if (sum == 0 && firstValueCount == 0) {
            firstValue = x;
            sum = 1;
            firstValueCount = 1;
            return;
        }
        if (sum == firstValueCount) {
            if (x == firstValue) {
                firstValueCount++;
                sum++;
                return;
            }
            xMin = firstValue;
            xMax = x;
            if (xMin > xMax) {
                xMax = firstValue;
                xMin = x;
            }
            deltaX = (xMax - xMin) / (nBins-1);
            xMin -= 0.5*deltaX;
            xMax += 0.5*deltaX;
            // resetting will clobber firstValueCount, so save it here
            long myFirstValueCount = firstValueCount;
            double myFirstValue = firstValue;
            reset();
            for (int i=0; i<myFirstValueCount; i++) {
                super.addValue(myFirstValue);
            }
            super.addValue(x);
            return;
        }
        // infinity will mess everything up, so just ignore it here
        if (!Double.isInfinite(x) && (x < xMin || x > xMax)) {
            collapseData(x);
        }
        super.addValue(x);
    }
    
    public void setXRange(DoubleRange xRange) {
        super.setXRange(xRange);
        // flag to prevent us from auto-scaling
        firstValueCount = -1;
    }
    
    /**
     * Collapses existing data and expands x range sufficiently so that the 
     * given x value is within the range of the histogram.  The collapsing is
     * done by taking the data from N bins and placing them into a single bin
     * over the whole range.  The resulting histogram is identical to the 
     * histogram that would have been collected if the resulting bin width had
     * been used the whole time.
     */
    protected void collapseData(double newX) {
        if (newX > xMin && newX < xMax) {
            // why were we called?
            return;
        }
        if (newX < xMin) {
            double newDeltaX = (xMax - newX) / nBins;
            int collapseFactor = (int)Math.ceil(newDeltaX/deltaX);
            long[] newCounts = new long[nBins];
            for (int i=0; i<nBins; i++) {
                newCounts[nBins-1-(i/collapseFactor)] += counts[nBins-1-i];
            }
            counts = newCounts;
            xMin -= (xMax - xMin) * (collapseFactor - 1);
            deltaX *= collapseFactor;
        }
        else {
            double newDeltaX = (newX - xMin) / nBins;
            int collapseFactor = (int)Math.ceil(newDeltaX/deltaX);
            long[] newCounts = new long[nBins];
            for (int i=0; i<nBins; i++) {
                newCounts[(i/collapseFactor)] += counts[i];
            }
            counts = newCounts;
            xMax += (xMax - xMin) * (collapseFactor - 1);
            deltaX *= collapseFactor;
        }
        
        for(int i=0; i<nBins; i++) {
            xValues[i] = xMin + (i+0.5)*deltaX;
        }
    }
    
    public double[] getHistogram() {
        if (sum > 0 && firstValueCount == sum) {
            // we have data, but not enough to set up the x range
            System.err.println("insufficient data to set up HistogramCollapsing, come back later");
        }
        return super.getHistogram();
    }
}
