/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.IMolecule;

/**
 * Interface for an iterator that can be targeted at one or more atoms.
 * Setting a target causes the iterator to produce only those of its iterates 
 * that contain the target atom(s).  A typical use would be to set an atom
 * as a target for a pair iterator, so the iterator would form only those pairs
 * that contain the targeted atom.  The argument can be null, indicating no 
 * target.
 */
public interface MoleculesetIteratorTargetable extends MoleculesetIterator {
	public void setTarget(IMolecule targetAtom);
}
