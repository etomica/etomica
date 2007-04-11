/* History
 * Created on Aug 4, 2004
 */
package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.atom.AtomSet;

/**
 * Interface for classes that loop over a set of atoms. Permits
 * iteration via a hasNext()-next() while loop (iterator returns
 * atoms to client) or via a call to allAtoms(AtomsetActive) (client gives
 * action to iterator).
 */

public interface AtomsetIterator {
    
    /**
     * Indicates whether the iterator has another atom.  
     * Once the iterator expires, this remains false, and will 
     * return true again only upon a call to reset.  No other methods
     * have the effect of making hasNext return true if it is
     * presently false.
     */
    public boolean hasNext();
    
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
     * Returns the next iterate without advancing the iterator.
     */
    public AtomSet peek();
    
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
