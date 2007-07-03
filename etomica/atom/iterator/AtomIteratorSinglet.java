package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.atom.AtomSet;
import etomica.atom.AtomSetSinglet;
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
    public AtomIteratorSinglet() {
        hasNext = false;
        atomSetSinglet = new AtomSetSinglet();
    }
    
    /**
     * Constructs iterator specifying that it return the given atom.  Call
     * to reset() must be performed before beginning iteration.
     * @param a The atom that will be returned by this iterator upon reset.
     */
    public AtomIteratorSinglet(IAtom a) {
        this();
        setAtom(a);
    }
        
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
    public int size() {return atom != null ? 1 : 0;}

	public void allAtoms(AtomsetAction action) {
		if (atom != null) {
            atomSetSinglet.atom = atom;
            action.actionPerformed(atomSetSinglet);
        }
	}
        
    /**
     * Sets iterator to a state where hasNext() returns false.
     */
    public void unset() {hasNext = false;}
    
    /**
     * Resets iterator to a state where hasNext is true.
     */
    public void reset() {
        hasNext = atom != null; 
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
        atomSetSinglet.atom = nextAtom();
        if (atomSetSinglet.atom == null) return null;
        return atomSetSinglet;
    }
    
    public final int nBody() {return 1;}
    
    private static final long serialVersionUID = 1L;
    private boolean hasNext = false;
    private IAtom atom;
    protected final AtomSetSinglet atomSetSinglet;
}
