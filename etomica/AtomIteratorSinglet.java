package etomica;

import etomica.action.AtomsetAction;

/**
 * Iterator that expires after returning a single atom, which is
 * specified by a call to the setAtom method, or via the constructor.
 * Subsequent calls to reset() and next() will return the specified atom,
 * until another is specified via setAtom.
 *
 * @author David Kofke
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
     * If atom is null, hasNext will remain false on reset.
     */
    public void setAtom(Atom a) {
    	atom[0] = a;
    	unset();
    }
    
    /**
     * returns 1.
     */
    public int size() {return 1;}

	public void allAtoms(AtomsetAction action) {
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
        hasNext = (atom[0] != null); 
    }
    
    /**
     * Returns the iterator's atom and unsets iterator.
     */
    public Atom nextAtom() {
    	if (!hasNext) return null;
    	hasNext = false;
    	return atom[0];
    }
    
    public Atom[] next() {
    	if (!hasNext) return null;
    	hasNext = false;
    	return atom;
    }
    
    /**
     * Returns the atom last specified via setAtom.  Does
     * not advance iterator.
     */
    public Atom[] peek() {
    	return hasNext ? null : atom;
    }
    
    public final int nBody() {return 1;}
    
    private boolean hasNext = false;
    private final Atom[] atom = new Atom[1];

}//end of AtomIteratorSinglet
