/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.api.IPotentialMolecular;
import etomica.molecule.IMoleculeList;

/**
 * Class defining a particular action to be performed on a set of atoms using an
 * arbitrary potential.  Examples of such actions are summing the energy, 
 * computing forces on atoms, determining collision times, etc.
 * Concrete subclasses define these actions through implementation of the 
 * doCalculation(IAtomSet, IPotential) method, which should perform the
 * defined calculation on the atoms using the given potential.
 *
 * @see PotentialMaster
 * @see PotentialGroup
 */
public interface PotentialCalculationMolecular extends PotentialCalculation {
 	
	/**
	 * Method giving the specific calculation performed by this class.
	 * @param atoms IAtomSet the atom sets for which the calculation is performed.
	 * @param potential The potential used to apply the action defined by this class.
	 */
	public void doCalculation(IMoleculeList atoms, IPotentialMolecular potential);
	
}
