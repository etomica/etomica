package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;

/**
 * Iterator that expires after returning a single atom, which is
 * specified by a call to the setAtom method, or via the constructor.
 * Subsequent calls to reset() and next() will return the specified atom,
 * until another is specified via setAtom.
 *
 * @author David Kofke
 */
public final class AtomIteratorSinglet implements AtomIteratorAtomDependent, java.io.Serializable {
    
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
    public AtomIteratorSinglet(IAtom a) {setAtom(a);}
        
    /**
     * Defines atom returned by iterator and leaves iterator unset.
     * Call to reset() must be performed before beginning iteration.
     * If atom is null, hasNext will remain false on reset.
     */
    public void setAtom(IAtom a) {
    	atom = a;
    	unset();
    }
    
    /**
     * @return the atom given by this iterator as its single iterate
     */
    public IAtom getAtom() {
        return atom;
    }
    
    /**
     * returns 1 if atom is not null, 0 if atom is null.
     */
    public int size() {return (atom == null) ? 0 : 1;}

	public void allAtoms(AtomsetAction action) {
		if(atom != null) action.actionPerformed(atom);
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
    public IAtom nextAtom() {
    	if (!hasNext) return null;
    	hasNext = false;
    	return atom;
    }
    
    public AtomSet next() {
        return nextAtom();
    }
    
    /**
     * Returns the atom last specified via setAtom.  Does
     * not advance iterator.
     */
    public AtomSet peek() {
    	return hasNext ? atom : null;
    }
    
    public final int nBody() {return 1;}
    
    private static final long serialVersionUID = 1L;
    private boolean hasNext = false;
    private IAtom atom;
}
