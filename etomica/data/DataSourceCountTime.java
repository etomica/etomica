package etomica.data;

import etomica.integrator.IntegratorMD;
import etomica.units.Time;

/**
 * Data source that keeps track of the elapsed simulation time of an MD
 * integrator. More precisely, sums the integrator's interval value times the
 * time-step each time the integrator fires an INTERVAL event. A START event
 * from the integrator will reset the elapsedTime.
 */
public class DataSourceCountTime extends DataSourceScalar {

    /**
	 * Sets up data source with no integrator specified.  Requires
	 * call to addIntegrator or setIntegrator before use.
	 */
	public DataSourceCountTime() {
		super("Simulation Time", Time.DIMENSION);
	}
    
    public DataSourceCountTime(IntegratorMD integrator) {
        this();
        setIntegrator(integrator);
    }

    public void setIntegrator(IntegratorMD newIntegrator) {
        integrator = newIntegrator;
    }
    
	/**
	 * Returns the simulation time elapsed by the integrator tracked
	 * by this class since the last reset. 
	 */
	public double getDataAsScalar() {
        return integrator.getCurrentTime();
	}
    
    private static final long serialVersionUID = 2L;
    protected IntegratorMD integrator;
}