package etomica;
 
/**
* Parent class for atom iterators.
* Atom iterators yield a sequence of atoms in successive calls to the next() method.
* When all atoms have been returned, hasNext() returns false.
* Iterators are often defined to progress "Up" or "Down" the set of atoms.
* "Up" and "Down" are arbitrary designations, except that the iterators guarantee
* that if atom 1 is "up list" of atom 2, then atom 2 is "down list" of atom 1.
* "Up" and "Down" relation between any atoms may change during the course of the 
* simulation, but at any instant the order is consistent and reproducible.
* "Neighbor" iterators yield only atoms that are considered to be "neighbors" of
* a specified atom.  The definition of "neighbor" depends on the iterator.  "Up neighbors"
* are those neighbors uplist of the atom; likewise with "Down neighbors".
*
* @see IteratorFactory
* @author David Kofke
*/
public abstract class AtomIterator implements java.io.Serializable {        

    protected boolean hasNext;
    
    /**
     * Iterator is constructed not ready for iteration.  Must call a reset
     * method before use.  hasNext returns false until then.
     */
    public AtomIterator() {
        hasNext = false;
    }
    
    /**
     * @return true if the iterator will return another atom with a subsequent 
     * call to next(), false otherwise.
     */
    public boolean hasNext() {return hasNext;}

    public void reset(IteratorDirective id) {
/*        IteratorDirective.Bounds bounds = id.bounds();
        if(bounds == IteratorDirective.ALL)             reset();
        else if(bounds == IteratorDirective.FIRST)      reset(id.firstAtom());
        else if(bounds == IteratorDirective.FIRST_LAST) reset(id.firstAtom(), id.lastAtom());
        else hasNext = false;*/
        switch(id.atomCount()) {
            case 0:  reset(); 
                     break;
            case 1:  reset(id.atom1()); 
                     break;
            case 2:  reset(id.atom1(), id.atom2()); 
                     break;
            default: hasNext = false; 
                     break;
        }
    }
    
    /**
     * Indicates if the given atom is among the iterates returned by this iterator.
     *
     * @return <code>true</code> if atom is one of the iterates.
     */
    public abstract boolean contains(Atom atom);

    /**
     * @return the next atom in the list
     */
    public abstract Atom next();
    
    /**
     * Resets the iterator, so that it is ready to go through its list again, beginning
     * with its natural first iterate (the identity of which depends on the iterator).
     *
     * @return the atom that will be returned with the first call to next() 
     */
    public abstract Atom reset();

    /**
     * Resets the iterator in reference to the given atom.  For some iterators, this means
     * that if the given atom it among its iterates, it will be the next one returned; otherwise
     * the iterator will give hasNext false.  For neighbor-type iterators, the next atom 
     * will be an appropriate neighbor of the given atom, if one exists.
     * Iterations proceed until the natural termination of the sequence.
     * 
     * @param first  the nominal first atom returned by the iterator
     * @return       the atom that will be returned with the next subsequent call to next()
     */
    public abstract Atom reset(Atom first);

    /**
     * Resets iterator in reference to the given atoms.  Initiation of the iteration is
     * as given by a call to reset(first), and termination of the iteration is done in reference
     * to the second given atom.  If this atom is among the iterates of this iterator, iteration
     * will terminate when this atom is reached and returned by next(); if it is not one of
     * the iterates, iteration terminates when the natural end of the iteration is reached.
     *
     * @param first the nominal first atom returned by the iterator
     * @param last  the nominal last atom returned by the iterator
     * @return      the atom that will be returned with the first call to next()
     */
    public abstract Atom reset(Atom first, Atom last);
    /**
     * Performs the given Action on each atom in the list in sequence.
     * 
     * @param act
     * @see Atom.Action
     */
    public abstract void allAtoms(AtomAction act);
    
    /**
     * A placeholder iterator that contains no atoms and always returns false for hasNext.
     */
    public static final AtomIterator NULL = new Null();
    private static final class Null extends AtomIterator {
        private Null() {super();}
        public Atom reset() {return null;}
        public Atom reset(Atom a) {return null;}
        public Atom reset(Atom a1, Atom a2) {return null;}
        public boolean contains(Atom a) {return false;}
        public void allAtoms(AtomAction act) {}
        public Atom next() {return null;}
    }//end of Null    

}//end of AtomIterator
    
