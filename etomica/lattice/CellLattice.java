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
     * @param D the dimension of the lattice
     * @param siteFactory makes the sites of the lattice 
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
    private final int[] idx;//a work array

    /**
     * Extends the SimpleLattice neighbor iterator to provide methods that
     * specify the neighbor range in terms of a distance, rather than an index range. 
     */
    public static class NeighborIterator extends SimpleLattice.NeighborIterator {
        public NeighborIterator(int D) {
            super(D);
            idx = new int[D];
        }
        
        /**
         * Defines neighbors to be all cells that are within the given distance
         * from a central cell.  Specifically, the integer range in direction
         * j is 1 + (int)(neighborDistance/cellSize[j]). 
         */
        public void setRange(double neighborDistance) {
            for(int i=0; i<D; i++) {
                idx[i] = 1+(int)(neighborDistance/((CellLattice)lattice).cellSize[i]);
            }
            super.setRange(idx);
        }

        /**
         * Defines neighbors to be all cells that are within the given distance
         * from a central cell, possible with a different distance specified in 
         * each direction.  Specifically, the integer range in direction
         * j is 1 + (int)(neighborDistance/cellSize[j]). 
         */
        public void setRange(double[] neighborDistance) {
            if(neighborDistance.length != D) throw new IllegalArgumentException("Dimension of neighborDistance array inconsistent with dimension of lattice");
            for(int i=0; i<D; i++) {
                idx[i] = 1+(int)(neighborDistance[i]/((CellLattice)lattice).cellSize[i]);
            }
            super.setRange(idx);
        }
        
        private final int[] idx;//a work array
    }//end of NeighborIterator
}
