/*
 * History
 * Created on Aug 6, 2004 by kofke
 */
package etomica.util;

/**
 * @author kofke
 *
 * Interface for a class that can tabulate a one-dimensional histogram 
 * of data values.
 */
public interface Histogram {

	public void addValue(double x);
	
	public double[] getHistogram();
	
	public double[] xValues();
	
	public int getCount();
	
	public void setNBins(int n);
	public int getNBins();
	
	public void setAutoScale(boolean b);
	public boolean isAutoScale();
	
	public void setXRange(DoubleRange range);
	public DoubleRange getXRange();
	
    /**
     * resets all histogram values and counts to zero
     */
	public void reset();
	
	public interface Factory {
		public Histogram makeHistogram();
		public Histogram makeHistogram(int n);
		public Histogram makeHistogram(int n, DoubleRange xRange);
	}
	
}
