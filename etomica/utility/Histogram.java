package etomica.utility;

import etomica.DataSource;
import etomica.units.Dimension;
//import ptolemy.plot.Plot;

public class Histogram implements DataSource.X {
	private double deltaX;
	private int  sum;
	private int[] counts;
	private double[] histogram;
    private double xvalues[];
//    private Plot plot;  //Histogram.Plot inner class
    private boolean autoScale;
    private boolean firstValue = true;
    private double xMin,xMinOld;
    private double xMax,xMaxOld ;
    private int nValues ;
    private String name, label, xLabel;
    private Dimension dimension = Dimension.NULL;
    private Dimension xDimension = Dimension.NULL;

    /**
     * Default constructor, making a 100-bin histogram.  Sets autoscale to true.
     */
    public Histogram() {this(100);}
    /**
     * Makes a new histogram instance having a number of bins given by the argument.
     * Sets autoscale to true.
     */
    public Histogram(int n) {
	    nValues = n;
	    counts = new int[n];
	    histogram = new double[n];
        xvalues = new double[n];
        autoScale = true;
    }    
    /**
     * Makes a new histogram instance, with n bins and limiting x values of
     * x0 and x1.  Sets autoscale to false.
     */
    public Histogram(int n, double x0, double x1) {
        this(n);
        autoScale = false;
	    xMin = x0;
	    xMax = x1;
	    reset();
	}
	
	public void setName(String s) {name = s;}
	public String toString() {return name;}
	
	public void setLabel(String s) {label = s;}
	public String getLabel() {return label;}
	
    public Dimension getDimension() {return dimension;}
    public void setDimension(Dimension dim) {dimension = dim;}
        
	public void setXLabel(String s) {xLabel = s;}
	public String getXLabel() {return xLabel;}
	
    public Dimension getXDimension() {return xDimension;}
    public void setXDimension(Dimension dim) {xDimension = dim;}
        
	public boolean isAutoScale() {return autoScale;}
	public void setAutoScale(boolean b) {autoScale = b;}
	
	public void reset() {                //resets all histogram values and counts to zero
	    sum = 0;
	    deltaX = (xMax-xMin)/(double)(nValues);
	    for(int i=0; i<nValues; i++) {
	        histogram[i] = 0.0; 
	        counts[i] = 0;
	        xvalues[i] = xMin + (i+0.5)*deltaX;
	    }
	    firstValue = autoScale;
	}
	
	public void addValues(double[] x) {  //updates histogram through addition of multiple new values
	    for(int j=0; j<x.length; j++) {addValue(x[j]);}
	}
	
/*    public void addValue(double x) {     //takes new value and updates histogram
	    int i = (int)Math.floor(((x-xMin)/deltaX));
	    i = (i >= nValues) ? nValues-1 : i;  // i = Max(i,nValues)
	    i = (i < 0) ? 0 : i;                 // i = Max(i,0);
	    counts[i]++;
	    sum++;
    }
 */   
    public void addValue(double x) {     //takes new value and updates histogram
        if(firstValue) {
            xMin = x - Math.abs(0.1*x);
            xMax = x + Math.abs(0.1*x);
            reset();
            firstValue = false;
        }
        int i;
        if(x < xMin) {
            xMinOld = xMin;
            if(autoScale) {xMin = x; redistribute();}
            i = 0;
        }
        else if(x > xMax) {
            xMaxOld = xMax;
            if(autoScale) {xMax = x; redistribute();}
            i = nValues-1;
        }
	    else {
	        i = (int)Math.floor(((x-xMin)/deltaX));
	        if(i == nValues){i--;}
	    }
	    counts[i]++;
	    sum++;
    }
  
    
    private void redistribute() {
        double[] newX = new double[nValues];
        double newDeltaX = (xMax - xMin)/(double)(nValues);
        for(int i=0; i<nValues; i++) {newX[i] = xMin + (i+0.5)*newDeltaX;}
        int[] newCounts = new int[nValues];
        int j = 0;
        newCounts[0] = 0;
         for(int i=0; i<nValues; i++) {
            if(j == nValues) break;
            if(newX[i]+0.5*newDeltaX < xvalues[j]-0.5*deltaX) continue;
            while(j < nValues && xvalues[j]+0.5*deltaX <= newX[i]+0.5*newDeltaX) {
                newCounts[i] += counts[j];
                j++;
            }
            if(j == nValues) break;
         //  if(xvalues[j] >= 0.0){
            int rem = (int)(counts[j]*((newX[i]+0.5*newDeltaX)-(xvalues[j]-0.5*deltaX))/deltaX);
            newCounts[i] += rem;
            if(i+1 < nValues) newCounts[i+1] = counts[j]-rem;
            j++;
          /* }
           else{// if(xvalues[j] < 0.0){
            int rem = (int)(counts[j]*(Math.abs(xvalues[j]-0.5*deltaX)-Math.abs(newX[i]+0.5*newDeltaX))/deltaX);
            newCounts[i] += rem;
            if(i+1 < nValues) newCounts[i+1] = counts[j]-rem;
            j++;
         }*/
        }
        deltaX = newDeltaX;
        counts = newCounts;
        xvalues = newX;
    }
    	    
    public double[] getHistogram() {      //returns an array representing the present histogram
        if(sum == 0) {return histogram;}
	    for(int i=0; i<nValues; i++) {
	        histogram[i] = (double)counts[i]/((double)sum*deltaX);
	    }
	    return histogram;
    }
    public double[] values(DataSource.ValueType dummy) {
        return getHistogram();
    }
    public double[] xValues() {return xvalues;}
    
    public etomica.DataSource.ValueType[] dataChoices() {return null;}
    
    public double[] getX() {return xvalues;}
    
/*    public Plot getPlot() {                //returns a Plot capable of displaying the histogram
        if(plot == null) plot = new Plot();
        return plot;
    }
    
    public class Plot extends ptolemy.plot.Plot {
        
        public Plot() {
            setTitle("Histogram");
            setYRange(0, 1);
            setXRange(xMin, xMax);
//            setPointsPersistence(nPoints);
            setImpulses(true); 
            setMarksStyle("none");
            update();
        }
        
        public void update() {  //redraws plot with current values of histogram
            clear(false);
            repaint();
            setXRange(xMin, xMax);
            double[] hist = getHistogram();
            double[] x = Histogram.this.getX();
            for(int i=0; i<hist.length; i++) {
//                System.out.println(x[i]+"  "+hist[i]);
                addPoint(0, x[i],hist[i],false);
            }
            repaint();
        }
    }
 */   
}
