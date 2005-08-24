package etomica.data;

import etomica.data.types.DataDouble;
import etomica.units.Dimension;
import etomica.util.NameMaker;

/**
 * Particular data source for which the data is a simple scalar of type double.
 */
 
public abstract class DataSourceScalar implements DataSource, java.io.Serializable {
    
    public DataSourceScalar(String label, Dimension dimension) {
        data = new DataDouble(label, dimension);
        setName(NameMaker.makeName(this.getClass()));
    }
    
    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }
    
    /**
     * Returns a single scalar value as the measurement for the given phase.
     * Subclasses define this method to specify the measurement they make.
     */
	public abstract double getDataAsScalar();
	
    
    /**
     * Causes the single getDataAsScalar(Phase) value to be computed and
     * returned for the given phase. In response to a getData() call,
     * MeterAbstract superclass will loop over all phases previously specified
     * via setPhase and collect these values into a vector and return them in
     * response to a getData() call.
     */
	public final DataDouble getDataDouble() {
		data.x = getDataAsScalar();
		return data;
	}
    
    public final Data getData() {
        return getDataDouble();
    }
	
    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }
	public final DataDouble data;
    private String name;
}
