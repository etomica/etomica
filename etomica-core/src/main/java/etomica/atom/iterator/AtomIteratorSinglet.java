/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.AtomSetSinglet;

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
    
    public IAtomList next() {
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
