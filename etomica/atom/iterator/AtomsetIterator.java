package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.atom.AtomSet;

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
	 * Returns the next AtomSet iterate, or null if hasNext() is false.
	 */
    public AtomSet next();
    
    /**
     * Performs given action over all the iterates of this iterator in its
     * current state.  The state of hasNext() is not relevant to this method,
     * and no call to reset() is needed beforehand. The set of atoms encountered 
     * in the allAtoms call are the same set that would be returned by 
     * looping using hasNext/next.
     */
    public void allAtoms(AtomsetAction action);

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
