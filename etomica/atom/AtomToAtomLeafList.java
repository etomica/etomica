package etomica.atom;

import etomica.api.IAtomLeaf;
import etomica.api.IAtomList;


/**
 * Interface for class that determines an AtomArrayList given an atom.
 */
public interface AtomToAtomLeafList {

    /**
     * Returns the ArrayList that this instance associates with the given atom.
     * Should return null if the atom is null.
     */
    public IAtomList getAtomList(IAtomLeaf atom);
}