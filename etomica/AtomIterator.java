/* History
 * Created on Aug 4, 2004
 */
package etomica;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public interface AtomIterator extends AtomSetIterator {
    
    public boolean contains(Atom atom);
    
    /**
     * Resets the iterator to loop through its iterates again.
     * The specific behavior depends on the type of iterator.
     */
    public void reset();
    
    /**
     * Puts iterator in a state in which next() returns null.
     */
    public void unset();
    
    /**
     * Returns the next atom in the iteration sequence.  It returns
     * null to indicate the end of iteration.  All subsequent calls
     * will return null until reset() is called.
     */
    public Atom next();
    
    /**
     * Performs given actions over all the iterates of this iterator.  Iterates
     * are defined according to most recent call to a reset method, or to
     * a class-dependent default if reset was not previously called.  A call to reset
     * followed by this method should cause iteration over the same set of atoms
     * that would be returned by looping using hasNext/next.  Unlike the
     * hasNext/next iteration, allAtoms can be called successively without intervening
     * reset calls to loop over the same set of atoms repeatedly.
     */
    public void all(Atom basis, IteratorDirective id, AtomActive action);

    /**
     * The number of iterates returned by this iterator, if iterating after
     * a call to the no-argument reset().  Some iterators give a value that
     * differs from this definition.  Neighbor iterators, for example, return
     * the total number of atoms in the basis, not the number of neighbors.
     */
    public int size(); 

}
