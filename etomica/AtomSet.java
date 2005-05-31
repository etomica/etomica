package etomica;


/**
 * Interface for classes that represent a (typically small) collection
 * of atoms.  Implemented in particular by Atom and AtomPair.
 */

/*
 * History
 * Created on Feb 18, 2005 by kofke
 */
public interface AtomSet {

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
     * is null.  Subclasses should override method taking Object as argument, 
     * check for instanceof AtomSet and pass on to this method.
     */
    public boolean equals(AtomSet atoms);
    
    /**
     * A convenient instance of a zero-length atom set.
     */
    public static final AtomSet NULL = new AtomSet() {
        public Atom getAtom(int i) {
            throw new IllegalArgumentException("Cannot get an atom from a NULL AtomSet");
        }
        public int count() {
            return 0;
        }
        public boolean equals(Object object) {
            if(!(object instanceof AtomSet)) return false;
            return equals((AtomSet)object);
        }
        public boolean equals(AtomSet atoms) {
            return atoms.count() == 0;
        }
    };
    
}
