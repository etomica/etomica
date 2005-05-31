package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomPair;
import etomica.AtomPairIterator;
import etomica.AtomSet;
import etomica.action.AtomsetAction;

/**
 * Iterator that expires after returning a single atom pair, which is
 * specified by a call to the setPair methods, or via the constructor.
 * Subsequent calls to reset() and next() will return the specified atoms,
 * until another is specified via setPair.
 *
 * @author David Kofke
 */
 
public final class ApiSinglet implements AtomPairIterator {
    
    /**
     * Constructs iterator without defining atom.  No atoms will
     * be given by this iterator until a call to setAtom is performed.
     */
    public ApiSinglet() {
        pair = new AtomPair();
        hasNext = false;
    }
    
    /**
     * Constructs iterator specifying that it return the given atoms in a pair.  Call
     * to reset() must be performed before beginning iteration.
     */
    public ApiSinglet(Atom a0, Atom a1) {
        this();
        setPair(a0, a1);
    }

    /**
     * Constructs iterator specifying that it return the given pair.  Call
     * to reset() must be performed before beginning iteration.
     */
    public ApiSinglet(AtomPair pair) {
        setPair(pair);
    }
            
    /**
     * Defines atoms returned by iterator and leaves iterator unset.
     * Call to reset() must be performed before beginning iteration.
     */
    public void setPair(Atom a0, Atom a1) {
        if(pair == null) pair = new AtomPair();
    	pair.atom0 = a0;
        pair.atom1 = a1;
    	unset();
    }

    /**
     * Defines pair returned by iterator and leaves iterator unset.
     * The given instance (not a copy) is returned by the iterator.
     * If given null, hasNext will always return false, even after reset.
     */
    public void setPair(AtomPair pair) {
        this.pair = pair;
    }
    
    /**
     * @return the pair given by this iterator as its single iterate
     */
    public AtomPair getPair() {
        return pair;
    }
    
    /**
     * returns 1, the number of iterates that would be returned on reset.
     */
    public int size() {return (pair == null) ? 0 : 1;}

	public void allAtoms(AtomsetAction action) {
		if(pair!= null) action.actionPerformed(pair);
	}
        
    /**
     * Returns true if the given atom equals the atom passed to the last call to setAtom(Atom).
     */
    public boolean contains(AtomSet a) {
    	return (a != null && a.equals(pair));
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
        hasNext = (pair != null); 
    }
    
    /**
     * Returns the iterator's pair and unsets iterator.
     */
    public AtomPair nextPair() {
    	if (!hasNext) return null;
    	hasNext = false;
    	return pair;
    }
    
    public AtomSet next() {
        return nextPair();
    }
    
    /**
     * Returns the atom last specified via setAtom.  Does
     * not advance iterator.
     */
    public AtomSet peek() {
    	return hasNext ? pair : null;
    }
    
    /**
     * Returns 2, indicating that this is a pair iterator.
     */
    public final int nBody() {return 2;}
    
    private boolean hasNext = false;
    private AtomPair pair;

}//end of ApiSinglet
