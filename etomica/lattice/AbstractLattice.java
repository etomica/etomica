package etomica.lattice;

/**
 * Interface for a generic lattice, which is a collection of sites that can be
 * accessed individually via specification of a set of integers.  Any object can
 * play the role of a site. 
 */
public interface AbstractLattice {

    /**
     * Dimension of the lattice.  The value of D describes the number of
     * integers needed to specify a site (via the site method).
     */
    public int D();

    /**
     * Returns the site specified by the given index.  The number of integers
     * in the index array must be equal to D, the lattice dimension.  Returned
     * object may be built and configured when called, or it may be pre-built
     * when array is constructed or sized and returned with call.  In either case,
     * this method should return a unique instance when invoked with different
     * arguments; invoking with the same argument may or may not return the same
     * instance, depending on how concrete lattice is implemented.
     */
    public Object site(int[] index);
    
}