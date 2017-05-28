/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.histogram;

import etomica.math.DoubleRange;

/**
 * Interface for a class that can tabulate a one-dimensional histogram 
 * of data values.
 * 
 * @author kofke
 */
public interface Histogram {

    /**
     * Adds the given value to the Histogram
     */
	public void addValue(double x);

    /**
     * Returns the histogram as an array of doubles
     */
	public double[] getHistogram();
	
	/**
     * Returns the x values from the histogram as an array of doubles
	 */
    public double[] xValues();
	
    /**
     * Returns the number of times the Histogram's addValue method was called.
     */
    public long getCount();

    /**
     * Sets the number of bins used by the Histogram.  The existing histogram
     * may be dropped or redistributed to the new bins.
     */
	public void setNBins(int n);
    
    /**
     * Returns the number of bins in the histogram.
     */
	public int getNBins();

    /**
     * Sets the range of x values used in the histogram (min and max).  The 
     * existing histogram may be dropped or redistributed to the new bins.
     */
	public void setXRange(DoubleRange range);
    
    /**
     * Returns the range of x values (min and max) used in the histogram.
     */
	public DoubleRange getXRange();
	
    /**
     * resets all histogram values and counts to zero
     */
	public void reset();
}
