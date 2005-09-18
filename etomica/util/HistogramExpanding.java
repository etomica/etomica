package etomica.util;


/**
 * Histogram that expands the X range to accomodate any values outside
 * the current range maintaining the bin width at a fixed value.  The
 * number of bins and the X range can be set explicitly, but only if the
 * change would not cause old data to be dropped.  
 *
 *  @author Andrew Schultz
 */
public class HistogramExpanding implements Histogram, java.io.Serializable {

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
     * Constructs a new instance with the give deltaX bin width and
     * range of X values.
     */
    public HistogramExpanding(double deltaX, DoubleRange xRange) {
        this.deltaX = deltaX;
        setXRange(xRange);
    }
    /**
     * Makes a new histogram instance, with n bins covering the
     * given range of X values.
     */
    public HistogramExpanding(int n, DoubleRange xRange) {
        deltaX = (xRange.maximum() - xRange.minimum()) / n;
        setXRange(xRange);
    }
	
    public boolean isAutoScale() {return true;}
    public void setAutoScale(boolean b) {
        if (!b) throw new IllegalArgumentException("HistogramExpanding always auto scales");
    }
    
    public void reset() {
        sum = 0;
        for(int i=0; i<nBins; i++) {
            counts[i] = 0;
        }
    }

    public void calcXValues() {
        for(int i=0; i<nBins; i++) {
            xvalues[i] = xMin + (i+0.5)*deltaX;
        }
    }    

    public void addValue(double x) {     //takes new value and updates histogram
        int i;
        if(x < xMin) {
            setXRange(new DoubleRange(x,xMax));
            i = 0;
        }
        else if(x > xMax) {
            setXRange(new DoubleRange(xMin,x));
            i = nBins-1;
        }
        else {
            i = (int)Math.floor(((x-xMin)/deltaX));
            if(i == nBins){i--;}
        }
        counts[i]++;
        sum++;
    }
    
    public void setXRange(DoubleRange xRange) {
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
        xvalues = new double[nBins];
        calcXValues();
    }
    public DoubleRange getXRange() {return new DoubleRange(xMin, xMax);}
      
    public int getNBins() {return nBins;}

    /**
     * Sets the number of bins to n.  Adds or removes bins from the end of 
     * the X range.  For more flexibility, call setXRange directly.   
     */
    public void setNBins(int n) {
        if (n == nBins) return;
        // subtract a half to prevent roundup to n+1
        setXRange(new DoubleRange(xMin,xMin+deltaX*(n-0.5)));
    }

    public double[] getHistogram() {
        if(sum != 0) {
            for(int i=0; i<nBins; i++) {
                histogram[i] = counts[i]/(sum*deltaX);
            }
        }
        return histogram;
    }
    
    public int getCount() {return sum;}
 
    public double[] xValues() {return xvalues;}
    
    public static final Histogram.Factory FACTORY = new Factory(1.0);
    
    public static class Factory implements Histogram.Factory {
        public Factory(double deltaX) {
            this.deltaX = deltaX;
        }
        public Histogram makeHistogram() {return new HistogramExpanding(deltaX);}
        public Histogram makeHistogram(int n) {
            return new HistogramExpanding(deltaX, new DoubleRange(0.0,deltaX*n));
        }
        // drop n on the floor. 
        public Histogram makeHistogram(int n, DoubleRange xRange) {return new HistogramExpanding(deltaX, xRange);}
        private final double deltaX;
    }

    private final double deltaX;
    private int  sum;
    private int[] counts;
    private double[] histogram;
    private double xvalues[];
    private double xMin;
    private double xMax;
    private int nBins;
}
