package etomica.lattice;


/**
 * Interface for an iterator in which the number of indices
 * in each dimension may be specified.
 */

/*
 * History
 * Created on Jan 6, 2005 by kofke
 */
public interface IndexIteratorSizable extends IndexIterator {

    /**
     * Specifies maximum index in each dimension (for dimension
     * i, maximum index is size[i]-1).
     */
    public void setSize(int[] size);
}
