package etomica.data;

import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.units.Dimension;

/**
 * Particular data source for which the data is a simple scalar of type double.
 */
 
public abstract class DataSourceScalar implements DataSource, java.io.Serializable {
    
    public DataSourceScalar(String label, Dimension dimension) {
        data = new DataDouble();
        dataInfo = new DataInfoDouble(label, dimension);
        tag = new DataTag();
        dataInfo.addTag(tag);
    }
    
    public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
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
	
	protected final DataDouble data;
    protected final IDataInfo dataInfo;
    protected final DataTag tag;
}
