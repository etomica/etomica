package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomIterator;
import etomica.AtomSet;
import etomica.action.AtomsetAction;
import etomica.utility.Arrays;

/**
 * Iterates over all the atoms given across an array of iterators.  
 *
 * @author David Kofke
 */
 
public final class AtomIteratorCompound implements AtomIterator {
    
    /**
     * Construct iterator to loop over all iterates obtained from each iterator 
     * in the given array.
     */
    public AtomIteratorCompound(AtomIterator[] iterators) {
        iteratorSet = iterators;
        reset();
    }
    
    /**
     * Adds the given iterator to the set of iterators collected by this 
     * compound iterator.  No check is made the iterator was not already added;
     * if an iterator is added more than once, it will be iterated for each addition,
     * as if it were a different instance each time.
     */
    public void addIterator(AtomIterator iterator) {
        if(iterator == null) return;
        iteratorSet = (AtomIterator[])Arrays.addObject(iteratorSet, iterator);
        unset();
    }
    
    /**
     * Removes the given iterator from the set of iterators collected by
     * this compound iterator.  If iterator was not previously added (or is
     * null), no action is taken.
     */
    public void removeIterator(AtomIterator iterator) {
        iteratorSet = (AtomIterator[])Arrays.removeObject(iteratorSet, iterator);
    }
    
    public boolean hasNext() {return hasNext;}
    
    public void unset() {hasNext = false;}
    
    public int nBody() {return 1;}
    
    public int size() {
        if(iteratorSet == null) return 0;
        int count = 0;
        for(int i=0; i<iteratorSet.length; i++) {
            count += iteratorSet[i].size();
        }
        return count;
    }
    
    public boolean contains(AtomSet atoms) {
        if(iteratorSet == null) return false;
        for(int i=0; i<iteratorSet.length; i++) {
            if(iteratorSet[i].contains(atoms)) return true;
        }
        return false;
    }
        
    public void reset() {
        if(!hasIterator()) {
            hasNext = false; 
            return;
        }
        index = 0;
        Atom next = null;
        
        iteratorSet[index].reset();
        while(!iteratorSet[index].hasNext() && index+1 < iteratorSet.length) {
        	index++;
        	iteratorSet[index].reset();
        }

        hasNext = iteratorSet[index].hasNext();
    }
    
    public AtomSet peek() {
    	return hasNext ? iteratorSet[index].peek() : null;
    }
    
    public AtomSet next() {
        return nextAtom();
    }//end of next

    public Atom nextAtom() {
        if(!hasNext) return null;
        Atom atom = iteratorSet[index].nextAtom();
        while(!iteratorSet[index].hasNext()) {
            if(++index < iteratorSet.length) {
                iteratorSet[index].reset();
            }
            else {
                hasNext = false;
                break;
            }
        }
        return atom;
    }
    
    public void allAtoms(AtomsetAction action) {
        for(int i=0; i<iteratorSet.length; i++) {
            iteratorSet[i].allAtoms(action);
        }
    }

    private boolean hasIterator() {
        return (iteratorSet.length > 0);
    }
    
    private AtomIterator[] iteratorSet = new AtomIterator[0];
    private boolean hasNext;
    private int index;
    
}//end of AtomIteratorCompound