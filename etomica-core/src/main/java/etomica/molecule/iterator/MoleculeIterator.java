/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule.iterator;

import etomica.molecule.IMolecule;

/**
 * Interface for classes that loop over a set of atoms. Permits
 * iteration via a hasNext()-next() while loop (iterator returns
 * atoms to client) or via a call to allAtoms(AtomActive) (client gives
 * action to iterator).
 */

public interface MoleculeIterator extends MoleculesetIterator {
                    
	/**
	 * Returns the next atom in the iteration sequence, or
     * null if hasNext() is false.
	 */
    public IMolecule nextMolecule();
}
