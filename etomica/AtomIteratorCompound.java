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
     * Constructs iterator that loops over all the children of the given
     * atom group's children (i.e., the group's grandchildren).  Registers
     * as a listener to the group so that iterator is updated if children are
     * added to the given group.  This constructor is used to form molecule
     * iterators, which loop over the children (molecules) of all species
     * agents in a phase (which are the child atoms of the phase's species master).
     */      
    public AtomIteratorCompound(Atom group) {
        this(makeIteratorArray(group));
        basis = group;
        group.node.parentPhase().speciesMaster.addListener(this);//sets iterator to update if group has addition or removal of child atoms
    }
    /**
     * Construct iterator to loop over all iterates obtained from each iterator 
     * in the given array.
     */
    public AtomIteratorCompound(AtomIterator[] iterators) {
        iteratorSet = iterators;
//        reset();
    }
    
    private static AtomIterator[] makeIteratorArray(Atom group) {
        if(group == null || group.node.childAtomCount() == 0) return null;
        AtomIteratorSequential groupIterator = new AtomIteratorSequential(group);
        AtomIterator[] iterators = new AtomIterator[group.node.childAtomCount()];
        groupIterator.reset();
        int i = 0;
        while(groupIterator.hasNext()) {
            iterators[i++] = new AtomIteratorSequential(groupIterator.next());
        }
        return iterators;
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
    
    public void setBasis(Atom a) {
        if(!(a instanceof AtomGroup)) {
            iteratorSet = null;
            basis = null;
            return;
        }
        iteratorSet = makeIteratorArray((AtomGroup)a);
        if(a == basis) return;//don't do this before makeIteratorArray, because this method can be called to update the iterator array for the same basis
        if(basis != null) basis.node.parentPhase().speciesMaster.removeListener(this);
        a.node.parentPhase().speciesMaster.addListener(this);
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