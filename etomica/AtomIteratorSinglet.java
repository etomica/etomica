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
 
 /* History of changes
  * 8/4/02 (DAK) Modified reset(Atom) to set basis to given atom while putting iterator ready for iteration
  *              Change made while attempting to enable operation of PistonCylinder
  * 8/5/02 (DAK) Commented out modification of 8/4/02, restoring to previous version.
  */
public class AtomIteratorSinglet extends AtomIterator {
    
    private Atom atom;
    private boolean hasNext;
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

	public void all(AtomSet basis, IteratorDirective id, final AtomSetActive action) {
		 if(!(basis instanceof Atom && action instanceof AtomActive)) return;
		 all((Atom)basis, id, (AtomActive)action);
	}
    
	public void all(Atom basis, IteratorDirective id, final AtomActive action) {
		if(basis == null) return;
		if(id.atomCount()==0 || id.atom1().node.isDescendedFrom(basis)) action.actionPerformed(basis);
	} 
    
    /**
     * Returns true if the given atom is the atom passed to the last call to setAtom(Atom).
     */
    public boolean contains(Atom a) {return (a != null && a.node.isDescendedFrom(atom));}
    
    public boolean hasNext() {return hasNext;}
    
    public void unset() {hasNext = false;}
    
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
     * Resets iterator to return basis atom if it is descended from the given atom.
     */
    public Atom reset(Atom a) {
/*   //     atom = a;
        hasNext = (atom != null);
        return atom;
*/        if(atom == null) hasNext = false;
        else hasNext = contains(a);
        return hasNext ? atom : null;
    }
        
    public Atom next() {hasNext = false; return atom;}
    
    public void allAtoms(AtomAction act) {
        if(atom == null) return;
        act.actionPerformed(atom);
    }
}//end of AtomIteratorSinglet
        
