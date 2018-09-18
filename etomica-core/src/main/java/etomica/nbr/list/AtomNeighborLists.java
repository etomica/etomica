/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.list;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.AtomArrayList;

/**
 * Class used to maintain neighbor lists.  Holds lists of atoms that were 
 * elsewhere deemed to be neighbors of the another atom.  A separate 
 * AtomArrayList is kept for each potential (potential => AtomArrayList mapping
 * is the responsibility of the consumer). 
 */
public class AtomNeighborLists {

    private static final AtomArrayList[] EMPTY_ATOMLIST_ARRAY = new AtomArrayList[0];
    protected AtomArrayList[] upList, downList;
	
    /**
     * Constructs sequencer for the given atom.
     */
    public AtomNeighborLists() {
        upList = EMPTY_ATOMLIST_ARRAY;
        downList = EMPTY_ATOMLIST_ARRAY;
    }
    
    /**
     * Adds the given atom as a "down" neighbor interacting via the potential 
     * with the given index.  
     * @param a the new downlist neighbor atom
     * @param index the of the potential between the atoms
     */
    public void addUpNbr(IAtom a, int index) {
        upList[index].add(a);
    }

    /**
     * Adds the given atom as a "down" neighbor interacting via the potential 
     * with the given index.  
     * @param a the new downlist neighbor atom
     * @param index the of the potential between the atoms
     */
    public void addDownNbr(IAtom a, int index) {
        downList[index].add(a);
    }

    /**
     * Returns an array of uplist-neighbor-atom lists.  Each list in the
     * array corresponds to a specific potential. A zero-length list indicates
     * that no concrete potentials apply to the atom.
     */
    public IAtomList[] getUpList() {
        return upList;
    }
	
    /**
     * Returns an array of downlist-neighbor-atom lists.  Each list in the
     * array corresponds to a specific potential. A zero-length list indicates
     * that no concrete potentials apply to the atom.
     */
    public IAtomList[] getDownList() {
        return downList;
    }
    
    /**
     * Sets the number of up and down lists maintained by this instance.
     */
    protected void setCapacity(int newCapacity) {
        if (newCapacity == upList.length) {
            return;
        }
        upList = new AtomArrayList[newCapacity];
        downList = new AtomArrayList[newCapacity];
        for (int i=0; i<newCapacity; i++) {
            upList[i] = new AtomArrayList(1);
            downList[i] = new AtomArrayList(1);
        }
    }
	
    /**
     * Clears neighbor lists, removing all listed neighbor atoms.
     */
	public void clearNbrs() {
		int length = upList.length;
		for (int i=0; i<length; i++) {
			upList[i].clear();
			downList[i].clear();
		}
	}
}
