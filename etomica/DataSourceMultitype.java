/*
 * Created on Jul 25, 2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica;

/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public interface DataSourceMultitype extends DataSource {

    public double[] getData(DataType type);
    
    /**
     * The types of values available from this data source.  One of these
     * is expected as the argument to the getData method.
     */
    public DataType[] dataChoices();

}
