package etomica;

import etomica.units.Dimension;

/**
 * A consumer of data.  Data goes in and might (or might not) come out.
 */
public interface DataSink {
    
    /**
     * Gives data to DataSink for processing, display, or whatever it does.
     */
	public abstract void putData(double[] values);
    
    /**
     * Applies a label used to describe the data when presented, 
     * for example as a plot legend.
     */
    public void setLabel(String label);
    
    /**
     * Sets the physical dimensions (e.g., length) of the data.
     */
    public void setDimension(Dimension dimension);

}