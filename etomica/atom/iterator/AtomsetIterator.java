/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;


/**
 * Interface for classes that loop over a set of atoms. Permits
 * iteration via a next()!=null while loop (iterator returns
 * atoms to client) or via a call to allAtoms(AtomsetActive) (client gives
 * action to iterator).
 */

public interface AtomsetIterator {
    
    /**
     * Resets the iterator to loop through its iterates again.
     */
    public void reset();
    
	/**
	 * Puts iterator in a state in which hasNext() returns false.
	 */
    public void unset();
    
    /**
     * The number of iterates returned by this iterator in its current state.
     */
    public int size(); 

    /**
     * Returns the number of atoms given in each iterate, i.e., the
     * size of the atom array returned with each call to next().  
     * @return
     */
    public int nBody();
}
