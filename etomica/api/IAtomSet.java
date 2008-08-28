package etomica.api;

/**
 * Interface for a set of IAtoms.  The IAtomSet might contain 0, 1, 2 or many
 * IAtoms.
 * 
 * @author Andrew Schultz
 */
public interface IAtomSet {

    /**
     * Returns the i-th atom, with numbering beginning from 0. 
     * If i is greater than count-1, throws an IllegalArgumentException.
     */
    public IAtom getAtom(int i);

    /**
     * @return the number of atoms in the set
     */
    public int getAtomCount();
}