package etomica;

/**
 * Iterates over all the atoms given across an array of iterators.  
 *
 * @author David Kofke
 */
 
public final class AtomIteratorCompound extends AtomIterator {
    
    private AtomIterator[] iteratorSet;
    private boolean hasNext;
    private final IteratorDirective directive = new IteratorDirective();
    private int index;
    private Atom basis = null;
    
    /**
     * Construct iterator to loop over all iterates obtained from each iterator 
     * in the given array.
     */
    public AtomIteratorCompound(AtomIterator[] iterators) {
        iteratorSet = iterators;
        reset();
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
    public Atom getBasis() {return basis;}
    
    public int size() {
        if(iteratorSet == null) return 0;
        int count = 0;
        for(int i=0; i<iteratorSet.length; i++) {
            count += iteratorSet[i].size();
        }
        return count;
    }
    
    public boolean contains(Atom atom) {
        if(iteratorSet == null) return false;
        for(int i=0; i<iteratorSet.length; i++) {
            if(iteratorSet[i].contains(atom)) return true;
        }
        return false;
    }
        
    public Atom reset() {
        return reset(directive.clear());
    }

    public Atom reset(IteratorDirective id) {
        if(iteratorSet == null || iteratorSet.length == 0) {hasNext = false; return null;}
        directive.copy(id);
        index = 0;
        Atom next = null;
        
        do next = iteratorSet[index].reset(id);
        while(next == null && ++index < iteratorSet.length);

        hasNext = (next != null);
        return next;
    }
    
    
    public Atom next() {
        Atom next = iteratorSet[index].next();
        while(!iteratorSet[index].hasNext()) {
            if(++index < iteratorSet.length) {
                iteratorSet[index].reset(directive);
            }
            else {
                hasNext = false;
                break;
            }
        }
        return next;
    }//end of next
    
    public void allAtoms(AtomAction act) {
        for(int i=0; i<iteratorSet.length; i++) {
            iteratorSet[i].allAtoms(act);
        }
    }
    
}//end of AtomIteratorCompound