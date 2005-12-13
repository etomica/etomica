package etomica.lattice;

/**
 * Lattice or crystal without cubic symmetry, which can be characterized by three
 * lattice constants, a, b, &c.
 * 
 * @author nancycribbin
 *
 */
public interface NonCubicLattice extends SpaceLattice {
    
    /**
     * Sets the size of the unit cell, and the lengths of each edge.
     * @param the lengths of each edge, in double array form.
     */
    public void setLatticeConstants(double[] lengths);
    
    /**
     * Accesses the length of each edge of a unit cell.
     */
    public double getLatticeConstant();
    
    
}
