/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.integrator.IntegratorHard;
import etomica.units.dimensions.Time;

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
