/*
 * Created on Aug 5, 2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public interface AtomPairIterator extends AtomSetIterator {
	public abstract void all(Atom basis, IteratorDirective id, AtomPairActive action);
	
	public abstract void all(AtomPair basis, IteratorDirective id, AtomPairActive action);
    
    public abstract void setBasis(Atom a1, Atom a2);
    
    /**
     * Returns the number of pairs that would be given by this iterator
     * after a call to the no-argument reset() method.
     */
    public abstract int size();        
    
    public abstract boolean hasNext();
    
    public abstract void reset(IteratorDirective id);
    
    
   /**
     * Resets the iterator, so that it is ready to go through all of its pairs.
     */
    public abstract void reset();
        
        
    public abstract AtomPair next();

    /**
     * Performs the given action on all pairs returned by this iterator.
     */
    public abstract void allPairs(AtomPairAction act);
    
}
