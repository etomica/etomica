/*
 * Created on Aug 5, 2004
 *
 */
package etomica;

/**
 * Interface for iterator that gives AtomPair instances.
 */
public interface AtomPairIterator {
	
	/**
	 * Performs given action on all pairs given by this iterator
	 * in its current state.  Does not require reset() before invoking.
	 * Any hasNext/next iteration in progress is subject to disruption.
	 * Iterates encountered by action in this method are the same
	 * as those given by hasNext/next. 
	 */
	public abstract void allPairs(AtomPairActive action);
	
	/**
	 * Indicates whether the given atom pair will be returned by the
	 * iterator during its iteration. The order of the atoms in the pair
	 * is significant (this means that a value of true is returned only if
	 * one of the pairs returned by the iterator will have the same two 
	 * atoms in the same atom1/atom2 position as the given pair). Not
	 * dependent on state of hasNext.
	 */
	public boolean contains(AtomPair pair);
	
	/**
	 * Returns the next iterate without advancing the iterator.
	 * Returns null if hasNext is false.
	 */
	public AtomPair peek();
	
	/**
	 * Puts iterator in state such that hasNext is false.
	 */
	public void unset();
	
    /**
     * Returns the number of pairs that would be given by this iterator
     * after a call to the no-argument reset() method. Does not require
     * hasNext() to be currently true.  Any iteration that might be
     * in progress in is subject to disruption (method does not maintain
     * iteration state of iterator).
     */
    public abstract int size();        
    
    /**
     * Indicates whether iterator has another iterate.
     */
    public abstract boolean hasNext();
    
	/**
	 * Resets the iterator, so that it is ready to go through 
	 * all of its pairs.
	 */
    public abstract void reset();
    
    /**
     * Returns the next iterate.  If called when hasNext is null,
     * no guaranteed behavior is specified.
     */
    public abstract AtomPair next();
}
