/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.integrator.IntegratorMD;
import etomica.units.dimensions.Time;

/**
 * Data source that returns the elapsed simulation time of an MD integrator.
 */
public class DataSourceCountTime extends DataSourceScalar {

	public DataSourceCountTime(IntegratorMD integrator) {
		super("Simulation Time", Time.DIMENSION);
		this.integrator = integrator;
	}

	/**
	 * Returns the simulation time elapsed by the integrator tracked
	 * by this class since the last reset.
	 */
	public double getDataAsScalar() {
		return integrator.getCurrentTime();
	}

	protected final IntegratorMD integrator;
}
