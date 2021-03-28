/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.exception.MethodNotImplementedException;

/**
 * interface for spherical 2-body potentials
 */
public interface Potential2Spherical extends IPotentialAtomic {
	/**
	 * The pair energy u(r^2) with no truncation applied.
	 * @param r2 the square of the distance between the particles.
	 */
	default double u(double r2) {
		throw new MethodNotImplementedException();
	}

}
