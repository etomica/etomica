package etomica.lattice;

/**
 * Interface for an iterator in which the number of indices
 * in each dimension may be specified.
 */

public interface IndexIteratorSizable extends IndexIterator {

    /**
     * Specifies maximum index in each dimension (for dimension
     * i, maximum index is size[i]-1).
     */
    public void setSize(int[] size);
}
