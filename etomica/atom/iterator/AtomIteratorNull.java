package etomica.atom.iterator;

import java.io.Serializable;

import etomica.api.IAtomLeaf;
import etomica.api.IAtomList;

/**
 * Static iterator that returns no atoms.
 * @author kofke
 */
public final class AtomIteratorNull implements AtomIterator, Serializable {

    // prevent instantiation.  Consumers should use the INSTANCE field.
    private AtomIteratorNull() {}
    
    public IAtomList next() {return null;}

    public IAtomLeaf nextAtom() {return null;}

    public void reset() {}

    public int size() {return 0;}

    public void unset() {}

    public int nBody() {return 1;}
    
    private static final long serialVersionUID = 1L;
    public static final AtomIteratorNull INSTANCE = new AtomIteratorNull();
}