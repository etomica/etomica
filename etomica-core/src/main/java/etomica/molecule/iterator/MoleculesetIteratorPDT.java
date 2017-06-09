/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule.iterator;

import etomica.atom.iterator.AtomsetIteratorDirectable;

/**
 * This is an interface for iterators that are box dependent, directable, 
 * and targetable.  This interface defines no new methods, but collects the 
 * appropriate interfaces into a single interface.
 */
public interface MoleculesetIteratorPDT extends MoleculesetIteratorBoxDependent,
        AtomsetIteratorDirectable, MoleculesetIteratorTargetable {

}
