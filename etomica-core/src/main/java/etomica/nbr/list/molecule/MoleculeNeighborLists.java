/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.list.molecule;

import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeArrayList;

/**
 * Class used to maintain neighbor lists.  Holds lists of molecules that were 
 * elsewhere deemed to be neighbors of the another molecule.  A separate 
 * MoleculeArrayList is kept for each potential (potential => MoleculeArrayList mapping
 * is the responsibility of the consumer). 
 * 
 * @author Tai Boon Tan
 *
 */
public class MoleculeNeighborLists implements java.io.Serializable {

    private static final long serialVersionUID = 2L;
    protected MoleculeArrayList[] upList, downList;
	
    /**
     * Constructs sequencer for the given molecule.
     */
    public MoleculeNeighborLists() {
        upList = new MoleculeArrayList[0];
        downList = new MoleculeArrayList[0];
    }
    
    /**
     * Adds the given atom as a "down" neighbor interacting via the potential 
     * with the given index.  
     * @param a the new downlist neighbor molecule
     * @param index the of the potential between the molecules
     */
    public void addUpNbr(IMolecule a, int index) {
        upList[index].add(a);
    }

    /**
     * Adds the given atom as a "down" neighbor interacting via the potential 
     * with the given index.  
     * @param a the new downlist neighbor molecule
     * @param index the of the potential between the molecules
     */
    public void addDownNbr(IMolecule a, int index) {
        downList[index].add(a);
    }

    /**
     * Returns an array of uplist-neighbor-molecule lists.  Each list in the
     * array corresponds to a specific potential. A zero-length list indicates
     * that no concrete potentials apply to the atom.
     */
    public IMoleculeList[] getUpList() {
        return upList;
    }
	
    /**
     * Returns an array of downlist-neighbor-molecule lists.  Each list in the
     * array corresponds to a specific potential. A zero-length list indicates
     * that no concrete potentials apply to the atom.
     */
    public IMoleculeList[] getDownList() {
        return downList;
    }
    
    /**
     * Sets the number of up and down lists maintained by this instance.
     */
    protected void setCapacity(int newCapacity) {
        if (newCapacity == upList.length) {
            return;
        }
        upList = new MoleculeArrayList[newCapacity];
        downList = new MoleculeArrayList[newCapacity];
        for (int i=0; i<newCapacity; i++) {
            upList[i] = new MoleculeArrayList(1);
            downList[i] = new MoleculeArrayList(1);
        }
    }
	
    /**
     * Clears neighbor lists, removing all listed neighbor molecules.
     */
	public void clearNbrs() {
		int length = upList.length;
		for (int i=0; i<length; i++) {
			upList[i].clear();
			downList[i].clear();
		}
	}
}
