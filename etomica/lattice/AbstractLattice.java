package etomica.lattice;

/**
 * Interface for a generic lattice, which is a collection of sites that can be
 * accessed individually via specification of a set of integers.  Any object can
 * play the role of a site. 
 * Lattice is retangular in the sense that the size in one dimension does not depend on
 * the index in another dimension (e.g., it cannot be triangular).  
 */
public interface AbstractLattice {

    /**
     * Dimension of the lattice.  The value of D describes the number of
     * integers needed to specify a site in the array.
     * @return
     */
    public int D();

    /**
     * Returns an array containing all the sites in the lattice.
     */
    public Object[] sites();

    /**
     * Returns the site specified by the given index.  The number of integers
     * in the index array must be equal to D, the lattice dimension.
     */
    public Object site(int[] index);
    
    /**
     * Return the size of the lattice, specified by the number of sites in
     * each dimension.
     */
    public int[] getSize();

    /**
     * Sets the size of the lattice, specified by the number of sites in each dimension.
     * @param dim
     */
    public void setSize(int[] dim);
}