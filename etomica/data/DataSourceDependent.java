/*
 * History
 * Created on Jul 25, 2004 by kofke
 */
package etomica.data;

import etomica.DataSource;

/**
 * Interface for an object that depends on one or more data sources.
 */
public interface DataSourceDependent {
	
	public DataSource[] getDataSource();
	
	public void setDataSource(DataSource[] source);
	
	public etomica.utility.IntegerRange dataSourceCountRange();
}
