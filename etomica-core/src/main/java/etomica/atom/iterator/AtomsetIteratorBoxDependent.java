/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.box.Box;

/**
 * Interface indicating that an iterator can determine appropriate
 * atoms for iteration given an arbitrary box.  
 */
public interface AtomsetIteratorBoxDependent extends AtomLeafsetIterator {

    /**
     * Sets the Box to pull iterates from
     * @throws NullPointerException if the Box is null
     */
	public void setBox(Box box);
	
}
