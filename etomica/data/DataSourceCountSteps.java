package etomica.data;

import etomica.integrator.Integrator;
import etomica.units.Quantity;

/**
 * Data source that fronts the Integrator's step count as a piece of Data.
 */
public class DataSourceCountSteps extends DataSourceScalar implements  
        java.io.Serializable {

    /**
	 * Sets up data source to count integrator steps.
	 */
	public DataSourceCountSteps() {
        super("Integrator steps", Quantity.DIMENSION);
	}

    public DataSourceCountSteps(Integrator integrator) {
        this();
        setIntegrator(integrator);
    }
    
    public void setIntegrator(Integrator newIntegrator) {
        integrator = newIntegrator;
    }
    
    public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }

	/**
	 * Returns the number of steps performed by the integrator
	 */
	public double getDataAsScalar() {
		return integrator.getStepCount();
	}
    
    private static final long serialVersionUID = 2L;
    protected Integrator integrator;
}
