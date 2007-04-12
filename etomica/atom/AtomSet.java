package etomica.atom;


/**
 * Interface for classes that represent a (typically small) collection
 * of atoms.  Implemented in particular by Atom and AtomPair.
 */
public interface AtomSet extends java.io.Serializable {

    /**
     * Returns the i-th atom, with numbering beginning from 0. 
     * If i is greater than count-1, throws an IllegalArgumentException.
     */
    public IAtom getAtom(int i);
    
    /**
     * @return the number of atoms in the set
     */
    public int count();
}
