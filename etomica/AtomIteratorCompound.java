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
    private boolean setAsNeighbor = false;
    private int index;
    Atom basis = null;
    
    /**
     * @param iterator the iterator of the parent atoms.
     */      
    public AtomIteratorCompound(AtomGroup group) {
        this(makeIteratorArray(group));
    }
    public AtomIteratorCompound(AtomIterator[] iterators) {
        iteratorSet = iterators;
        reset();
    }
    
    private static AtomIterator[] makeIteratorArray(AtomGroup group) {
        if(group == null || group.childAtomCount() == 0) return null;
        AtomIteratorSequential groupIterator = new AtomIteratorSequential(group);
        AtomIterator[] iterators = new AtomIterator[group.childAtomCount()];
        groupIterator.reset();
        int i = 0;
        while(groupIterator.hasNext()) {
            iterators[i++] = new AtomIteratorSequential(groupIterator.next());
        }
        return iterators;
    }
    
    public boolean hasNext() {return hasNext;}
    
    public void setBasis(Atom a) {
        if(!(a instanceof AtomGroup)) {
            iteratorSet = null;
            basis = null;
            return;
        }
        iteratorSet = makeIteratorArray((AtomGroup)a);
        basis = a;
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
    
    /**
     * This method performs no action for this iterator.  Iterator
     * is fixed to always be NOT a neighbor iterator.
     */
    public void setAsNeighbor(boolean b) {
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
    
    //not implemented
    public void allAtoms(AtomAction act) {
    }
}//end of AtomIteratorCompound