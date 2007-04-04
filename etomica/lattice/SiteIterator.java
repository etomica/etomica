/*
 * History
 * Created on Dec 28, 2004 by kofke
 */
package etomica.lattice;


/**
 * Interface for objects that iterate over some or all of the
 * sites of a lattice.
 */
public interface SiteIterator {

    /**
     * Indicates whether the iterator has another site to return.
     */
    public boolean hasNext();

    
    /**
     * Returns the index of the next iterate while advancing iterator. A call to
     * peek() following by nextIndex() is the best way to obtain both the next
     * site and its index, if both are desired at once. Take care to copy values
     * in another array if planned for use after next iteration--values in
     * returned array may be changed when iterator is advanced. No meaningful
     * value is returned if hasNext is false.
     */
    public int[] nextIndex();
        
    /**
     * Returns the next iterate.
     */
    public Object next();

    /**
     * Returns the lattice index of the next iterate without advancing iterator.
     * No meaningful value is returned if hasNext is false.  Take care to copy
     * values in another array if planned for use with next iterate--values in
     * returned array may be changed when iterator is advanced.
     */
//    public int[] peekIndex();
    
    /**
     * Returns the next iterate without advancing the iterator.
     */
    public Object peek();

    /**
     * Puts the iterator in a state ready to begin iteration.
     */
    public void reset();

    /**
     * Returns the number of iterates that this iterator would give
     * if reset.
     */
    public int size();

    /**
     * Puts the iterator in a state such that hasNext is false.
     *
     */
    public void unset();
    
    /**
     * Sets the lattice having the sites to be given by this iterator.
     */
    public void setLattice(FiniteLattice lattice);
}
