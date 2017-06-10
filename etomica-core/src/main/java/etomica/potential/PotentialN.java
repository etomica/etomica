/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Space;

/**
 * @author kofke
 *
 * General potential that depends on positions of all N molecules, or is
 * otherwise not naturally expressed as a single-, pair-, etc-body potential.
 */

public abstract class PotentialN extends Potential {

	/**
	 * Constructor for PotentialN.
	 * @param space
	 */
	public PotentialN(Space space){
		super(Integer.MAX_VALUE, space);
	}
}
