package etomica.atom.iterator;

import etomica.AtomSet;
import etomica.AtomsetIterator;
import etomica.Debug;
import etomica.action.AtomsetAction;

/**
 * Iterator that expires after returning a single atom set, which is
 * specified by a call to the setAtom method, or via the constructor.
 * Subsequent calls to reset() and next() will return the specified atom,
 * until another is specified via setAtom.
 *
 * @author David Kofke
 */
 
 public class AtomsetIteratorSinglet implements AtomsetIterator, java.io.Serializable {
    
    /**
     * Constructs iterator without defining atom.  No atoms will
     * be given by this iterator until a call to setAtom is performed.
     */
    public AtomsetIteratorSinglet(int nBody) {
        this.nBody = nBody;
        hasNext = false;
    }
    
    /**
     * Constructs iterator specifying that it return the given atom.  Call
     * to reset() must be performed before beginning iteration.
     * @param a The atom that will be returned by this iterator upon reset.
     */
    public AtomsetIteratorSinglet(AtomSet a) {
    	this(a.count());
    	setAtom(a);
    }
        
    /**
     * Defines atom returned by iterator and leaves iterator unset.
     * Call to reset() must be performed before beginning iteration.
     */
    public void setAtom(AtomSet a) {
        if (Debug.ON && a != null && a.count() != nBody) throw new IllegalArgumentException("Wrong AtomSet count");
    	atom = a;
    	unset();
    }
    
    /**
     * returns 1 if atom has been specified, zero otherwise.
     */
    public int size() {return atom != null ? 0 : 1;}

	public void allAtoms(AtomsetAction action) {
		action.actionPerformed(atom);
	}
    
    /**
     * Returns true if the given atom equals the atom passed to the last call to setAtom(Atom).
     */
    public boolean contains(AtomSet a) {
    	return !(atom == null || a == null) && atom.equals(a);
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
        hasNext = true; 
    }
    
    /**
     * Returns the iterator's atom and unsets iterator.
     */
    public AtomSet next() {
    	if (!hasNext) return null;
    	hasNext = false;
    	return atom;
    }
    
    /**
     * Returns the atom last specified via setAtom.  Does
     * not advance iterator.
     */
    public AtomSet peek() {
    	return hasNext ? null : atom;
    }
    
    public final int nBody() {return nBody;}
    
    private final int nBody;
    private boolean hasNext = false;
    private AtomSet atom;

}//end of AtomsetIteratorSinglet
        
