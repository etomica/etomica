package etomica.util;

import java.io.Serializable;


/* History
 * 09/08/02 (DAK) added set/get methods for xMin, xMax, nBins
 * 08/04/04 (DAK,AJS,NRC) deleted DataSource.X methods; de-implemented DataSource.X.  Dimension-related material removed
 */

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
    private boolean firstValue;

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
        firstValue = true;
    }

	public void reset() {                //resets all histogram values and counts to zero
	    super.reset();
	    firstValue = true;
	}
	
    public void addValue(double x) {     //takes new value and updates histogram
        if(firstValue) {
            xMin = x - Math.abs(0.1*x);
            xMax = x + Math.abs(0.1*x);
            deltaX = (xMax-xMin)/nBins;
            reset();
            firstValue = false;
        }
        if (x < xMin || x > xMax) {
            collapseData(x);
        }
        super.addValue(x);
    }
    
    public void setXRange(DoubleRange xRange) {
        super.setXRange(xRange);
        firstValue = false;
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
            int[] newCounts = new int[nBins];
            for (int i=0; i<nBins; i++) {
                newCounts[nBins-1-(i/collapseFactor)] += counts[nBins-1-i];
            }
            xMin -= (xMax - xMin) * (collapseFactor - 1);
            deltaX *= collapseFactor;
        }
        else {
            double newDeltaX = (newX - xMin) / nBins;
            int collapseFactor = (int)Math.ceil(newDeltaX/deltaX);
            int[] newCounts = new int[nBins];
            for (int i=0; i<nBins; i++) {
                newCounts[(i/collapseFactor)] += counts[i];
            }
            xMax += (xMax - xMin) * (collapseFactor - 1);
            deltaX *= collapseFactor;
        }
        
        for(int i=0; i<nBins; i++) {
            xValues[i] = xMin + (i+0.5)*deltaX;
        }
    }
    
    public static final Histogram.Factory FACTORY = new Factory(100);
    
    public static class Factory implements Histogram.Factory, Serializable {
        public Factory(int n) {
            nBins = n;
        }
        
		public Histogram makeHistogram() {
            return new HistogramCollapsing(nBins);
        }
        
        private final int nBins;
    }
}
