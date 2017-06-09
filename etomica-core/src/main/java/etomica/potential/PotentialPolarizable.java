/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IMoleculeList;

/**
 * Interface for a polarizable potential.  The interface includes only a method
 * that returns the polarization energy from the last call to energy(AtomSet).
 * 
 * @author Ken
 */
public interface PotentialPolarizable {

    /**
     * Returns the polarization energy from the last call to energy(AtomSet).
     * @return
     */
	public double getLastPolarizationEnergy();

	/**
	 * Returns the polarization contribution (non pairwise-additive) to the
	 * energy for the given AtomSet.
	 */
	public double getPolarizationEnergy(IMoleculeList atoms);
}
