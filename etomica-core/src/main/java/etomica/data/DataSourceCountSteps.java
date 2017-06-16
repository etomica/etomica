/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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

	/**
	 * Returns the number of steps performed by the integrator
	 */
	public double getDataAsScalar() {
		return integrator.getStepCount();
	}
    
    private static final long serialVersionUID = 2L;
    protected Integrator integrator;
}
