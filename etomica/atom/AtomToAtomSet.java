package etomica.atom;


/**
 * Interface for class that determines an AtomArrayList given an atom.
 */
public interface AtomToAtomSet {

    /**
     * Returns the ArrayList that this instance associates with the given atom.
     * Should return null if the atom is null.
     */
    public AtomSet getAtomSet(IAtom atom);
}