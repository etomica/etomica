package etomica.nbr.list;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;

/**
 * Class used to maintain neighbor lists.  Holds lists of atoms that were 
 * elsewhere deemed to be neighbors of the another atom.  A separate 
 * AtomArrayList is kept for each potential (potential => AtomArrayList mapping
 * is the responsibility of the consumer). 
 */
public class AtomNeighborLists implements java.io.Serializable {

    private static final long serialVersionUID = 2L;
    protected AtomArrayList[] upList, downList;
	
    /**
     * Constructs sequencer for the given atom.
     */
    public AtomNeighborLists() {
        upList = new AtomArrayList[0];
        downList = new AtomArrayList[0];
    }
    
    /**
     * Adds the given atom as a "down" neighbor interacting via the potential 
     * with the given index.  
     * @param a the new downlist neighbor atom
     * @param index the of the potential between the atoms
     */
    public void addUpNbr(Atom a, int index) {
        upList[index].add(a);
    }

    /**
     * Adds the given atom as a "down" neighbor interacting via the potential 
     * with the given index.  
     * @param a the new downlist neighbor atom
     * @param index the of the potential between the atoms
     */
    public void addDownNbr(Atom a, int index) {
        downList[index].add(a);
    }

    /**
     * Returns an array of uplist-neighbor-atom lists.  Each list in the
     * array corresponds to a specific potential. A zero-length list indicates
     * that no concrete potentials apply to the atom.
     */
    public AtomArrayList[] getUpList() {
        return upList;
    }
	
    /**
     * Returns an array of downlist-neighbor-atom lists.  Each list in the
     * array corresponds to a specific potential. A zero-length list indicates
     * that no concrete potentials apply to the atom.
     */
    public AtomArrayList[] getDownList() {
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
            upList[i] = new AtomArrayList();
            downList[i] = new AtomArrayList();
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
