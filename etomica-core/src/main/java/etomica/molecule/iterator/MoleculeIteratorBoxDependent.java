/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule.iterator;


/**
 * Interface indicating that an atom iterator can determine appropriate
 * atoms for iteration given an arbitrary box.  
 */
public interface MoleculeIteratorBoxDependent extends MoleculeIterator, MoleculesetIteratorBoxDependent {

}
