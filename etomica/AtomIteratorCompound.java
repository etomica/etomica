package etomica;

/**
 * Iterates over all the atoms given across an array of iterators.  
 *
 * @author David Kofke
 */
 
public final class AtomIteratorCompound implements AtomIterator {
    
    private AtomIterator[] iteratorSet;
    private boolean hasNext;
    private final IteratorDirective directive = new IteratorDirective();
    private int index;
    private Atom atoms[];
    
    /**
     * Construct iterator to loop over all iterates obtained from each iterator 
     * in the given array.
     */
    public AtomIteratorCompound(AtomIterator[] iterators) {
        iteratorSet = iterators;
        reset();
        atoms = new Atom[1];
    }

	public void all(Atom basis, IteratorDirective id, final AtomActive action) {
		if(basis == null || basis.node.isLeaf() || action == null) return;
		throw new RuntimeException("Method all not implemented in AtomIteratorCompound");
	}
    
    /**
     * Adds the given iterator to the set of iterators collected by this 
     * compound iterator.
     */
    public void addIterator(AtomIterator iterator) {
        if(iterator == null) return;
        int currentLength = (iteratorSet != null) ? iteratorSet.length : 0;
        AtomIterator[] newSet = new AtomIterator[currentLength+1];
        for(int i=0; i<currentLength; i++) newSet[i] = iteratorSet[i];
        newSet[currentLength] = iterator;
        iteratorSet = newSet;
    }
    
    public boolean hasNext() {return hasNext;}
    
    public void unset() {hasNext = false;}
    
    //try to eliminate this method
    public void setBasis(Atom a) {
        throw new RuntimeException("AtomIteratorCompound.setBasis not defined");
    }
    public Atom getBasis() {return null;}
    
    public int nBody() {return 1;}
    
    public int size() {
        if(iteratorSet == null) return 0;
        int count = 0;
        for(int i=0; i<iteratorSet.length; i++) {
            count += iteratorSet[i].size();
        }
        return count;
    }
    
    public boolean contains(Atom[] atoms) {
        if(iteratorSet == null) return false;
        for(int i=0; i<iteratorSet.length; i++) {
            if(iteratorSet[i].contains(atoms)) return true;
        }
        return false;
    }
        
    public void reset() {
        if(iteratorSet == null || iteratorSet.length == 0) {hasNext = false; return;}
        index = 0;
        Atom next = null;
        
        iteratorSet[index].reset();
        while(!iteratorSet[index].hasNext() && index+1 < iteratorSet.length) {
        	index++;
        	iteratorSet[index].reset();
        }

        hasNext = iteratorSet[index].hasNext();
    }
    
    public Atom[] peek() {
    	return iteratorSet[index].peek();
    }
    
    public Atom[] next() {
        atoms[0] = iteratorSet[index].nextAtom();
        while(!iteratorSet[index].hasNext()) {
            if(++index < iteratorSet.length) {
                iteratorSet[index].reset();
            }
            else {
                hasNext = false;
                break;
            }
        }
        return atoms;
    }//end of next

    public Atom nextAtom() {
    	return next()[0];
    }
    
    public void allAtoms(AtomsetActive act) {
        for(int i=0; i<iteratorSet.length; i++) {
            iteratorSet[i].allAtoms(act);
        }
    }
    
}//end of AtomIteratorCompound