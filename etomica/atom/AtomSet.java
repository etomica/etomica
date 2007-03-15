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
    public Atom getAtom(int i);
    
    /**
     * @return the number of atoms in the set
     */
    public int count();
    
    /**
     * Element-by-element comparison of equality of this atom set with
     * another.  Order of atoms is relevant.  Returns false if given atom set
     * is null.  Returns true if all atoms in two sets are equal, even if not
     * comparing the same instance of AtomSet.
     * Subclasses should also override the method that takes Object as argument, 
     * check for instanceof AtomSet and pass on to this method.
     */
    public boolean equals(AtomSet atoms);
}
