package etomica;

/**
 * Iterates over all the atoms given across an array of iterators.  
 *
 * @author David Kofke
 */
 
public final class AtomIteratorCompound implements AtomIterator, PhaseListener {
    
    private AtomIterator[] iteratorSet;
    private boolean hasNext;
    private final IteratorDirective directive = new IteratorDirective();
    private boolean setAsNeighbor = false;
    private int index;
    Atom basis = null;
    
    /**
     * Construct iterator to loop over all iterates obtained from each iterator 
     * in the given array.
     */
    public AtomIteratorCompound(AtomIterator[] iterators) {
        iteratorSet = iterators;
//        reset();
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
    
    //try to eliminate this method
    public void setBasis(Atom a) {
        throw new RuntimeException("AtomIteratorCompound.setBasis not defined");
/*        if(!(a instanceof AtomGroup)) {
            iteratorSet = null;
            basis = null;
            return;
        }
        iteratorSet = makeIteratorArray((AtomGroup)a);
        if(a == basis) return;//don't do this before makeIteratorArray, because this method can be called to update the iterator array for the same basis
        if(basis != null) basis.node.parentPhase().speciesMaster.removeListener(this);
        a.node.parentPhase().speciesMaster.addListener(this);
        basis = a;*/
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
        throw new RuntimeException("Method allAtoms not implemented in AtomIteratorCompound");
    }
    
    /**
     * Updates set of iterators in the event that the basis group has an
     * addition or removal of one of its children.  This method is called
     * by the species master of the basis group when it fires an event
     * in its addNotify method.
     */
    public void actionPerformed(PhaseEvent evt) {
        if(basis == null) return;
        if(evt.type() == PhaseEvent.ATOM_ADDED || evt.type() == PhaseEvent.ATOM_REMOVED) {
            if(evt.atom().node.parentGroup() == basis) setBasis(basis);
        }
    }
    public void actionPerformed(SimulationEvent evt) {
        actionPerformed((PhaseEvent)evt);
    }
}//end of AtomIteratorCompound