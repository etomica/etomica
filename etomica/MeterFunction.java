package etomica;

/**
 * Meter for recording a function
 *
 * @author David Kofke
 */
public abstract class MeterFunction extends MeterArray implements DataSource {
    
    protected DataSource xDataSource;
    
    /**
     * Default constructor.  Sets the xDataSource to a DataSourceUniform
     */
    public MeterFunction() {
        this(new DataSourceUniform());
    }

    public MeterFunction(DataSource xDataSource) {
    	super();
    	this.xDataSource = xDataSource;
    }
	 
    /**
     * Returns the DataSource for the X values of this meter
     */
    public DataSource getXDataSource() {return xDataSource;}
        
}