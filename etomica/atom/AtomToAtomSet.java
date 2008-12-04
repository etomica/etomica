package etomica.atom;

import etomica.api.IAtom;
import etomica.api.IAtomList;


/**
 * Interface for class that determines an AtomArrayList given an atom.
 */
public interface AtomToAtomSet {

    /**
     * Returns the ArrayList that this instance associates with the given atom.
     * Should return null if the atom is null.
     */
    public IAtomList getAtomSet(IAtom atom);
}