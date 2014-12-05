/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.api.IBox;

/**
 * Interface indicating that an iterator can determine appropriate
 * atoms for iteration given an arbitrary box.  
 */
public interface MoleculesetIteratorBoxDependent extends MoleculesetIterator {

    /**
     * Sets the Box to pull iterates from
     * @throws a NullPointerException if the Box is null
     */
	public void setBox(IBox box);
	
}
