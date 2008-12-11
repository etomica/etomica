package etomica.api;

/**
 * Interface for a set of IAtoms.  The IAtomSet might contain 0, 1, 2 or many
 * IAtoms.
 * 
 * @author Andrew Schultz
 */
public interface IMoleculeList {

    /**
     * Returns the i-th atom, with numbering beginning from 0. 
     * If i is greater than count-1, throws an IllegalArgumentException.
     */
    public IMolecule getMolecule(int i);

    /**
     * @return the number of atoms in the set
     */
    public int getMoleculeCount();
}