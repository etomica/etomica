package etomica.api;

/**
 * Interface for a set of IAtoms.  The 
 * @author andrew
 *
 */
public interface IAtomSet {

    /**
     * Returns the i-th atom, with numbering beginning from 0. 
     * If i is greater than count-1, throws an IllegalArgumentException.
     */
    public abstract IAtom getAtom(int i);

    /**
     * @return the number of atoms in the set
     */
    public abstract int getAtomCount();
}