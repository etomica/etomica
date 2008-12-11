package etomica.atom;

import etomica.api.IMolecule;
import etomica.api.IMoleculeList;


/**
 * Interface for class that determines an AtomArrayList given an atom.
 */
public interface MoleculeToMoleculeList {

    /**
     * Returns the ArrayList that this instance associates with the given atom.
     * Should return null if the atom is null.
     */
    public IMoleculeList getMoleculeList(IMolecule atom);
}