package etomica.util;

import java.io.Serializable;


/* History
 * 09/08/02 (DAK) added set/get methods for xMin, xMax, nBins
 * 08/04/04 (DAK,AJS,NRC) deleted DataSource.X methods; de-implemented DataSource.X.  Dimension-related material removed
 */

public class HistogramSimple implements Histogram, java.io.Serializable {
	private double deltaX;
	private int  sum;
	private int[] counts;
	private double[] histogram;
    private double xvalues[];
    private boolean autoScale;
    private boolean firstValue = true;
    private double xMin;
    private double xMax;
 //   private DoubleRange xRange;
    private int nBins ;
    private String name;

    /**
     * Default constructor, making a 100-bin histogram.  Sets autoscale to true.
     */
    public HistogramSimple() {this(100);}
    /**
     * Makes a new histogram instance having a number of bins given by the argument.
     * Sets autoscale to true.
     */
    public HistogramSimple(int n) {
	    nBins = n;
	    counts = new int[n];
	    histogram = new double[n];
        xvalues = new double[n];
        autoScale = true;
    }    
    /**
     * Makes a new histogram instance, with n bins and limiting x values of
     * x0 and x1.  Sets autoscale to false.
     */
    public HistogramSimple(int n, DoubleRange xRange) {
        this(n);
        autoScale = false;
	    setXRange(xRange);
	}
	
	public void setName(String s) {name = s;}
	public String toString() {return name;}
	
	public boolean isAutoScale() {return autoScale;}
	public void setAutoScale(boolean b) {autoScale = b;}
	
	public void reset() {                //resets all histogram values and counts to zero
	    sum = 0;
	    deltaX = (xMax-xMin)/(double)(nBins);
	    for(int i=0; i<nBins; i++) {
	        histogram[i] = 0.0; 
	        counts[i] = 0;
	        xvalues[i] = xMin + (i+0.5)*deltaX;
	    }
	    firstValue = true;
	}
	
	public void addValues(double[] x) {  //updates histogram through addition of multiple new values
	    for(int j=0; j<x.length; j++) {addValue(x[j]);}
	}
	
    public void addValue(double x) {     //takes new value and updates histogram
        if(firstValue && autoScale) {
            xMin = x - Math.abs(0.1*x);
            xMax = x + Math.abs(0.1*x);
            reset();
            firstValue = false;
        }
        int i;
        if(x < xMin) {
            if(autoScale) {xMin = x; redistribute();}
            i = 0;
        }
        else if(x > xMax) {
            if(autoScale) {xMax = x; redistribute();}
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
//        this.xRange = xRange;
        xMin = xRange.minimum();
        xMax = xRange.maximum();
        reset();
    }
    public DoubleRange getXRange() {return new DoubleRange(xMin, xMax);}
      
    public int getNBins() {return nBins;}
    public void setNBins(int n) {
        this.nBins = n;
        counts = new int[n];
        histogram = new double[n];
        xvalues = new double[n];
        reset();
    }
    
    private void redistribute() {
        double[] newX = new double[nBins];
        double newDeltaX = (xMax - xMin)/(double)(nBins);
        for(int i=0; i<nBins; i++) {newX[i] = xMin + (i+0.5)*newDeltaX;}
        int[] newCounts = new int[nBins];
        int j = 0;
        newCounts[0] = 0;
         for(int i=0; i<nBins; i++) {
            if(j == nBins) break;
            if(newX[i]+0.5*newDeltaX < xvalues[j]-0.5*deltaX) continue;
            while(j < nBins && xvalues[j]+0.5*deltaX <= newX[i]+0.5*newDeltaX) {
                newCounts[i] += counts[j];
                j++;
            }
            if(j == nBins) break;
         //  if(xvalues[j] >= 0.0){
            int rem = (int)(counts[j]*((newX[i]+0.5*newDeltaX)-(xvalues[j]-0.5*deltaX))/deltaX);
            newCounts[i] += rem;
            if(i+1 < nBins) newCounts[i+1] = counts[j]-rem;
            j++;
          /* }
           else{// if(xvalues[j] < 0.0){
            int rem = (int)(counts[j]*(Math.abs(xvalues[j]-0.5*deltaX)-Math.abs(newX[i]+0.5*newDeltaX))/deltaX);
            newCounts[i] += rem;
            if(i+1 < nBins) newCounts[i+1] = counts[j]-rem;
            j++;
         }*/
        }
        deltaX = newDeltaX;
        counts = newCounts;
        xvalues = newX;
    }
    	    
    public double[] getHistogram() {      //returns an array representing the present histogram
        if(sum != 0) {
		    for(int i=0; i<nBins; i++) {
		        histogram[i] = (double)counts[i]/((double)sum*deltaX);
		    }
        }
	    return histogram;
    }
    
    public int getCount() {return sum;}
 
    public double[] xValues() {return xvalues;}
    
    public static final Histogram.Factory FACTORY = new HistogramSimpleFactory();
    
    public static class HistogramSimpleFactory implements Histogram.Factory, Serializable {
		public Histogram makeHistogram() {return new HistogramSimple();}
		public Histogram makeHistogram(int n) {return new HistogramSimple(n);}
		public Histogram makeHistogram(int n, DoubleRange xRange) {return new HistogramSimple(n, xRange);}
    }
}
