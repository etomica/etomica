/*
 * Created on Aug 5, 2004
 *
 */
package etomica;

/**
 * @author andrew
 *
 */
public interface AtomPairIterator {
	public abstract void allPairs(AtomPairActive action);
	
	public boolean contains(AtomPair pair);
	
	public AtomPair peek();
	
	public void unset();
	
    /**
     * Returns the number of pairs that would be given by this iterator
     * after a call to the no-argument reset() method.
     */
    public abstract int size();        
    
    public abstract boolean hasNext();
    
   /**
     * Resets the iterator, so that it is ready to go through all of its pairs.
     */
    public abstract void reset();
        
    public abstract AtomPair next();
}
