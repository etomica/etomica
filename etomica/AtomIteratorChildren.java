package etomica;

/**
 * Iterates over the child atoms of a sequence of atoms.
 * Parent atoms themselves are not returned in the sequence.
 *
 * @author David Kofke
 */
 
public final class AtomIteratorChildren implements AtomIterator {
    
    private AtomIterator parentIterator;
    private boolean hasNext;
    private final IteratorDirective directive = new IteratorDirective();
    private boolean setAsNeighbor = false;
    private AtomIterator currentIterator;
    
    public AtomIteratorChildren(AtomIterator iterator) {
        parentIterator = iterator;
    }
    
    public boolean hasNext() {return hasNext;}
    
    public boolean contains(Atom atom) {
        parentIterator.reset();
        while(parentIterator.hasNext()) {
            if(((AtomGroup)parentIterator.next()).childIterator.contains(atom)) return true;
        }
        return false;
    }
    
    public Atom reset(IteratorDirective id) {
        directive.copy(id);
        parentIterator.reset(id);
        return reset();
    }
    
    public void setAsNeighbor(boolean b) {
        parentIterator.setAsNeighbor(b);
        setAsNeighbor = b;
    }
    
    public Atom reset() {
        Atom atom = null;
        parentIterator.setAsNeighbor(setAsNeighbor);
        parentIterator.reset();
        while(parentIterator.hasNext()) {
            currentIterator = ((AtomGroup)parentIterator.next()).childIterator;
            currentIterator.setAsNeighbor(setAsNeighbor);
            atom = currentIterator.reset(directive);
            if(currentIterator.hasNext()) {hasNext = true; break;}
        }
        return atom;
    }
    
    public Atom next() {
        Atom next = currentIterator.next();
        while(!currentIterator.hasNext()) {
            if(parentIterator.hasNext()) {
                currentIterator = ((AtomGroup)parentIterator.next()).childIterator;
                currentIterator.setAsNeighbor(setAsNeighbor);
                currentIterator.reset(directive);
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
}//end of AtomIteratorChildren