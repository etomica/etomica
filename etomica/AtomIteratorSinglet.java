package etomica;

/**
 * Iterator that expires after returning a single atom, which is
 * specified by a call to the setAtom method, or via the constructor.
 * reset() sets the iterator to return the atom.
 * reset(Atom) sets as follows:
 *   if isAsNeighbor is false, the iterator will return its atom only if
 *   the atom given to reset matches it.
 *   if isAsNeighbor is true, the iterator will return its atom only if
 *   the atom is appropriately up or down list from the given atom, in correspondence
 *   with the current value of the iterator's direction field.
 *
 * @author David Kofke
 */
public class AtomIteratorSinglet implements AtomIterator {
    
    private Atom atom;
    private boolean hasNext, isAsNeighbor;
    private IteratorDirective.Direction direction;
    
    public AtomIteratorSinglet() {hasNext = false;}
    public AtomIteratorSinglet(Atom a) {setAtom(a);}
        
    /**
     * Mutator method for this iterator's atom.
     */
    public void setAtom(Atom a) {atom = a; hasNext = false;}
    /**
     * Accessor method for this iterator's atom.
     */
    public Atom getAtom() {return atom;}
    
    public void setBasis(Atom a) {setAtom(a);}
    public Atom getBasis() {return atom;}
    
    public int size() {return (atom != null) ? 1 : 0;}
    
    /**
     * Returns true if the given atom is the atom passed to the last call to setAtom(Atom).
     */
    public boolean contains(Atom a) {return (a != null && a.isDescendedFrom(atom));}
    
    public boolean hasNext() {return hasNext;}
    
    public void setAsNeighbor(boolean b) {isAsNeighbor = b;}
    
    public Atom reset(IteratorDirective id) {
        direction = id.direction();
        switch(id.atomCount()) {
            case 0:  return reset(); 
            case 1:  return reset(id.atom1()); 
  //          case 2:  return reset(id.atom1(), id.atom2()); 
            default: hasNext = false; 
            return null;
        }
    }
    
    /**
     * Resets iterator to return the iterator's atom.  Ignores any specifications
     * of direction or isAsNeighbor.
     */
    public Atom reset() {
        hasNext = (atom != null); 
        return atom;
    }
    
    /**
     * Sets the given atom as the one returned by the iterator.
     */
    public Atom reset(Atom a) {
        if(atom == null) hasNext = false;
        else if(isAsNeighbor) {
            if(direction == IteratorDirective.UP) hasNext = a.preceeds(atom);
            else if(direction == IteratorDirective.DOWN) hasNext = atom.preceeds(a);
            else if(direction == IteratorDirective.BOTH) hasNext = true;
            else /*direction == NEITHER*/ hasNext = false;
        }
        else hasNext = contains(a);  // isAsNeighbor == false
        return hasNext ? atom : null;
    }
        
    public Atom next() {hasNext = false; return atom;}
    
    public void allAtoms(AtomAction act) {
        if(atom == null) return;
        act.actionPerformed(atom);
    }
}//end of AtomIteratorSinglet
        
