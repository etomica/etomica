/*
 * History
 * Created on Dec 17, 2004 by kofke
 */
package etomica.lattice;

import etomica.Space;


/**
 * A lattice of arbitrarily-sized rectangular cells, such that a point in space can be 
 * associated with one of the cells.
 */
public class CellLattice extends SimpleLattice {

    /**
     * @param D
     * @param siteFactory
     */
    public CellLattice(int D, SiteFactory siteFactory) {
        super(D, siteFactory);
        cellSize = new double[D];
        idx = new int[D];
    }
    
    /**
     * Returns the cell in which the given point lies.  An ArrayIndexOutOfBounds exception
     * is thrown if the point lies outside the bounds of the lattice (the boundary
     * in each direction j is given by cellSize[j]*dimensions[j]).  The dimension of
     * the vector is assumed to be consistent with the dimension D of the lattice (i.e.,
     * r.D() == this.D()) but this is not checked.
     */
    public Object site(Space.Vector r) {
        for(int i=0; i<idx.length; i++) {
            idx[i] = (int)(r.x(i)/cellSize[i]);
//            double x = r.x(i)/cellSize[i];
//            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
//            while(idx[i] >= dimensions[i]) idx[i] -= dimensions[i];
//            while(idx[i] < 0)              idx[i] += dimensions[i];
        }
        return site(idx);

    }

    /**
     * Specifies the spatial dimensions of all cells in the lattice, such that a cell in the
     * lattice will have a size in each direction as given by the corresponding value
     * in the argument.  An IllegalArgumentException is thrown if the length of the
     * given array is not equal to the dimension D() of the lattice.
     */
    public void setCellSize(double[] newSize) {
        if(newSize.length != this.cellSize.length) throw new IllegalArgumentException("Dimension of cell size array inconsistent with dimensions of cells");
        System.arraycopy(newSize, 0, this.cellSize, 0, newSize.length);
    }
    
    private final double[] cellSize;
    private final int[] idx;
}
