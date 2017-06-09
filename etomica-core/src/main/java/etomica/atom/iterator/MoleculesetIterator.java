/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.IMoleculeList;

/**
 * Interface for classes that loop over a set of atoms. Permits
 * iteration via a next()!=null while loop (iterator returns
 * atoms to client) or via a call to allAtoms(AtomsetActive) (client gives
 * action to iterator).
 */

public interface MoleculesetIterator extends AtomsetIterator {
    
	/**
	 * Returns the next AtomSet iterate, or null if hasNext() is false.
	 */
    public IMoleculeList next();
}
