/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.IAtom;

/**
 * Interface for an atom iterator that can be altered by setting
 * of an atom.  A neighbor iterator, for example. 
 */
public interface AtomIteratorAtomDependent extends AtomIterator {

	public void setAtom(IAtom atom);
}
