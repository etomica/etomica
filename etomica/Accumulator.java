/*
 * History
 * Created on Jul 20, 2004
 */
package etomica;

import etomica.units.Dimension;

/**
 * Parent class for all accumulators, which are responsible for gathering
 * and processing data during a simulation.
 * @author kofke
 */
public abstract class Accumulator implements DataSource, DataSink {

    /**
     * Number of points processed by accumulator. 
     */
    protected int nData;
    
    private String label = "Accumulator";
    private Dimension dataSourceDimension;
    
	/**
	 * No-argument constructor
	 */
	public Accumulator() {
		this(Dimension.UNDEFINED);
	}
	
	public Accumulator(Dimension dataSourceDimension) {
		this.dataSourceDimension = dataSourceDimension;
	}
	/**
	 * Accessor method for the meter's label
	 */
	public String getLabel() {return label;}
	/**
	 * Accessor method for the meter's label
	 */
	public void setLabel(String s) {label = s;}


	public abstract void add(double[] values);
	
	public abstract double[] getData();
	
	public abstract void reset();
	
	public Dimension getDimension() {return dataSourceDimension;}
	
	public void setDataSourceDimension(Dimension dataSourceDimension) {
		this.dataSourceDimension = dataSourceDimension;
	}
}
