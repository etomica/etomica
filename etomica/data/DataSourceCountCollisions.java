package etomica.data;

import etomica.integrator.IntegratorHard;
import etomica.units.Time;

/**
 * This is a data source to count the number of collisions processed by a
 * hard-potential integrator.
 */
public class DataSourceCountCollisions extends DataSourceScalar {

    /**
     * Sets up data source with no integrator specified.  Requires
     * call to addIntegrator or setIntegrator before use.
     */
    public DataSourceCountCollisions() {
        super("Collision Count", Time.DIMENSION);
    }

    public DataSourceCountCollisions(IntegratorHard integrator) {
        this();
        setIntegrator(integrator);
    }
    
    public void setIntegrator(IntegratorHard newIntegrator) {
        integrator = newIntegrator;
    }
    
    /**
     * Returns the simulation time elapsed by the integrator tracked
     * by this class since the last reset. 
     */
    public double getDataAsScalar() {
        return integrator.getCollisionCount();
    }
    
    private static final long serialVersionUID = 2L;
    protected IntegratorHard integrator;
}