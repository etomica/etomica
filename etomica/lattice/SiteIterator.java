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
     * Returns the next iterate.
     */
    public Object next();

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
    public void setLattice(AbstractLattice lattice);
}
