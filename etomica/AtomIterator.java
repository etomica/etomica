package etomica;

/**
* Interface for atom iterators,  yield a sequence of atoms in successive calls 
* to the next() method.
* When all atoms have been returned, hasNext() returns false.
* Iterators are often defined to progress "Up" or "Down" the set of atoms.
* "Up" and "Down" are arbitrary designations, except that the iterators guarantee
* that if atom 1 is "up list" of atom 2, then atom 2 is "down list" of atom 1.
* "Up" and "Down" relation between any atoms may change during the course of the 
* simulation, but at any instant the order is consistent and reproducible.
*
* @see IteratorFactory
* @author David Kofke
*/

public interface AtomIterator {
    
    public boolean hasNext();
    
    public boolean contains(Atom atom);
    
    public Atom reset(IteratorDirective id);
    
    public void setAsNeighbor(boolean b);
    
    public Atom reset();
    
    public Atom next();
    
    public void allAtoms(AtomAction act);
    

    /**
     * A placeholder iterator that contains no atoms and always returns false for hasNext.
     */
    public static final AtomIterator NULL = new Null();
    static final class Null implements AtomIterator {
        public boolean hasNext() {return false;}
        public boolean contains(Atom atom) {return false;}
        public Atom reset(IteratorDirective id) {return null;}
        public void setAsNeighbor(boolean b) {}
        public Atom reset() {return null;}
        public Atom next() {return null;}
        public void allAtoms(AtomAction act) {}
    }//end of Null    

    /**
     * An enumerated type used to specify if an atom iterator is to skip
     * the first atom of its natural sequence, or to include it.  Skipping
     * it would be appropriate if, for example, the iterator is providing the inner
     * loop of an atom pair iterator.  The Initiation field of an iterator derived
     * from AtomIteratorAbstract is final and cannot be changed after construction.
     */
/*    public static final class Initiation extends Constants.TypedConstant {
            
        private Initiation(String label) {super(label);}
        public static final Initiation[] CHOICES = new Initiation[] {
            new Initiation("Skip first"),
            new Initiation("Include first"),
        };
        
        public final Constants.TypedConstant[] choices() {return CHOICES;}
    }//end of Initiation
    public static final Initiation SKIP_FIRST = Initiation.CHOICES[0];
    public static final Initiation INCLUDE_FIRST = Initiation.CHOICES[1];
*/
}//end of AtomIterator