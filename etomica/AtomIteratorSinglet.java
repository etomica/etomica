package etomica;

/**
 * Iterator that expires after returning a single atom, which is
 * specified by a call to the reset(Atom) method.
 * A second atom may be specified using the reset(Atom, Atom) method. 
 * Then if the iterator is set to isAsNeighbor, the second atom
 * is the only one returned; otherwise
 * the first atom is the only one returned.
 *
 * @author David Kofke
 */
public class AtomIteratorSinglet implements AtomIterator {
    
    private Atom atom1, atom2, next;
    private boolean hasNext, isAsNeighbor;
    
    public AtomIteratorSinglet() {hasNext = false;}
    public AtomIteratorSinglet(Atom a) {reset(a);}
        
    /**
     * Returns true if the given atom is the atom passed to the last call to reset(Atom),
     * or is one of the atoms passed to the last call to reset(Atom, Atom).
     */
    public boolean contains(Atom a) {return (a == atom1 || a == atom2) && a != null;}
    
    public boolean hasNext() {return hasNext;}
    
    public void setAsNeighbor(boolean b) {isAsNeighbor = b;}
    
    public Atom reset(IteratorDirective id) {
        switch(id.atomCount()) {
            case 0:  return reset(); 
            case 1:  return reset(id.atom1()); 
            case 2:  return reset(id.atom1(), id.atom2()); 
            default: hasNext = false; 
            return null;
        }
    }
    
    /**
     * Resets iterator to return atom specified by previous call to reset(Atom), or
     * to the second atom of the last call to reset(Atom, Atom) if SKIP_FIRST.
     */
    public Atom reset() {
        next = isAsNeighbor ? atom1 : atom2;
        hasNext = (next != null); 
        return next;
    }
    
    /**
     * Sets the given atom as the one returned by the iterator if INCLUDE_FIRST.
     */
    public Atom reset(Atom a) {
        atom1 = a;
        atom2 = null;
        return reset();
    }
    
    /**
     * 
     */
    public Atom reset(Atom first, Atom last) {
        atom1 = first;
        atom2 = last;
        return reset();
    }
    
    public Atom next() {hasNext = false; return next;}
    
    public void allAtoms(AtomAction act) {
        reset();
        act.actionPerformed(next);
    }
}//end of AtomIteratorSinglet
        
