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
     * in the index array must be equal to D, the lattice dimension.  No specification
     * is made here regarding the identity of the instance returned by this method.
     * Thus repeated calls may return different object instances, or the same instance,
     * regardless of the argument.  The concrete implementations of this class should
     * provide this specification. 
     */
    public Object site(int[] index);
    
}