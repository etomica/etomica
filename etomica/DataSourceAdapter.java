/*
 * History
 * Created on Oct 21, 2004 by kofke
 */
package etomica;

import etomica.units.Dimension;

/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public abstract class DataSourceAdapter implements DataSource {

	private String label;
	private final Dimension dimension;
	
	public DataSourceAdapter(Dimension dimension) {
		this.dimension = dimension;
	}
	
	/* (non-Javadoc)
	 * @see etomica.DataSource#getLabel()
	 */
	public String getLabel() {
		return label;
	}
	
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

}
