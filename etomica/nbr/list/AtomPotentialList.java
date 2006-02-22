/*
 * History
 * Created on Sep 20, 2004 by kofke
 */
package etomica.nbr.list;

import java.io.Serializable;

import etomica.util.Arrays;

/**
 * Class used to maintain list of whether each 1-body potential that
 * applies to an AtomType is currently in effect or not for an Atom.
 */
public class AtomPotentialList implements Serializable {

    protected boolean[] isInteractingList;
	
    /**
     * Constructs sequencer for the given atom.
     */
    public AtomPotentialList() {
        isInteractingList = new boolean[0];
    }
    
    /**
     * Sets whether the atom is interacting or not with the potential
     * corresponding to the given index.
     */
    public void setIsInteracting(boolean flag, int index) {
        ensureCapacity(index);
        if (!(index < isInteractingList.length)) {
            isInteractingList = Arrays.resizeArray(isInteractingList,index+1);
        }
        isInteractingList[index] = flag;
    }

    /**
     * Returns booleans indicating whether the Atom interacts with each 
     * potential or not.
     */
    public boolean[] getInteractingList() {
        return isInteractingList;
    }
	
    protected void ensureCapacity(int index) {
        if (index > isInteractingList.length-1) {
            isInteractingList = Arrays.resizeArray(isInteractingList, index);
        }
    }
	
    /**
     * Should be called when removing a potential that applied to this atom
     */
	public void decrementInteractingList() {
		if (isInteractingList.length == 0) throw new RuntimeException("potential list empty in removePotential");
		isInteractingList = new boolean[isInteractingList.length-1];
	}
	
    /**
     * Clears neighbor lists, resetting all flags to false (not interacting).
     */
	public void clearNbrs() {
		final int length = isInteractingList.length;
		for (int i=0; i<length; i++) {
			isInteractingList[i] = false;
		}
	}
}
