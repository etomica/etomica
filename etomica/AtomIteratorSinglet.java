package etomica;

/**
 * Iterator that expires after returning a single atom, which is
 * specified by a call to the setAtom method, or via the constructor.
 * Subsequent calls to reset() and next() will return the specified atom,
 * until another is specified via setAtom.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 8/4/02 (DAK) Modified reset(Atom) to set basis to given atom while putting iterator ready for iteration
  *              Change made while attempting to enable operation of PistonCylinder
  * 8/5/02 (DAK) Commented out modification of 8/4/02, restoring to previous version.
  * 08/26/04 (DAK) revised with overhaul of iterators
  */
public final class AtomIteratorSinglet implements AtomIteratorAtomDependent {
    
    /**
     * Constructs iterator without defining atom.  No atoms will
     * be given by this iterator until a call to setAtom is performed.
     */
    public AtomIteratorSinglet() {hasNext = false;}
    
    /**
     * Constructs iterator specifying that it return the given atom.  Call
     * to reset() must be performed before beginning iteration.
     * @param a The atom that will be returned by this iterator upon reset.
     */
    public AtomIteratorSinglet(Atom a) {setAtom(a);}
        
    /**
     * Defines atom returned by iterator and leaves iterator unset.
     * Call to reset() must be performed before beginning iteration.
     */
    public void setAtom(Atom a) {
    	atom[0] = a;
    	unset();
    }
    
    /**
     * Returns 0 if atom has not been previously set, or has been set to null;
     * returns 1 otherwise.
     */
    public int size() {return (atom[0] != null) ? 1 : 0;}

	public void allAtoms(AtomsetActive action) {
		if(atom[0] != null) action.actionPerformed(atom);
	}
        
    /**
     * Returns true if the given atom equals the atom passed to the last call to setAtom(Atom).
     */
    public boolean contains(Atom[] a) {
    	return (a != null && a[0] != null && a[0].equals(atom));
    }
    
    /**
     * Returns true if the an atom has been set and a call to reset() has been
     * performed, without any subsequent calls to next().
     */
    public boolean hasNext() {return hasNext;}
    
    /**
     * Sets iterator to a state where hasNext() returns false.
     */
    public void unset() {hasNext = false;}
    
    /**
     * Resets iterator to a state where hasNext is true.
     */
    public void reset() {
        hasNext = (atom != null); 
    }
    
    /**
     * Returns the iterator's atom and unsets iterator.
     */
    public Atom nextAtom() {hasNext = false; return atom[0];}
    
    public Atom[] next() {hasNext = false; return atom;}
    
    /**
     * Returns the atom last specified via setAtom.  Does
     * not advance iterator.
     */
    public Atom[] peek() {
    	return atom;
    }
    
    public final int nBody() {return 1;}
    
    private boolean hasNext = false;
    private final Atom[] atom = new Atom[1];

}//end of AtomIteratorSinglet
        
