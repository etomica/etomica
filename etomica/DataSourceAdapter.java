/*
 * History
 * Created on Oct 21, 2004 by kofke
 */
package etomica;

import etomica.units.Dimension;

/**
 * Convenience class to permit easy implementation of DataSource
 * interface.  Defines methods to set and get labels, constructor
 * that defines the dimension, and a simple identity translator.
 * Subclasses need implement only the getData() method.  If a different
 * DataTranslator is appropriate, subclass should override the getTranslator method.
 */
public abstract class DataSourceAdapter implements DataSource {

	/**
	 * Constructs a data source that returns data of the given dimension.
	 * @param dimension physical dimension of the data from this source
	 */
	public DataSourceAdapter(Dimension dimension) {
		this.dimension = dimension;
	}
	
	/**
	 * Returns a string the describes the data given by this source.
	 */
	public String getLabel() {
		return label;
	}
	
	/**
	 * Sets the descriptive label for this source.
	 * @param label new value for the label
	 */
	public void setLabel(String label) {
		this.label = label;
	}

	/* (non-Javadoc)
	 * @see etomica.DataSource#getDimension()
	 */
	public Dimension getDimension() {
		return dimension;
	}

	/* (non-Javadoc)
	 * @see etomica.DataSource#getTranslator()
	 */
	public DataTranslator getTranslator() {
		return DataTranslator.IDENTITY;
	}
	
	private String label;
	protected Dimension dimension;
	
}
