/*
 * History
 * Created on Dec 17, 2004 by kofke
 */
package etomica.lattice;

import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.Default;

/**
 * A lattice of arbitrarily-sized rectangular cells, such that a point in space can be 
 * associated with one of the cells.
 */
public class CellLattice extends RectangularLattice {

    /**
     * @param dimensions
     *            the spatial dimensions of the lattice, such that the total
     *            space occupied by the lattice cells a rectangular box with
     *            each side of length dimensions.x(i)/size[i] for side i. The
     *            dimensions at any time are given by the array passed here, so
     *            any changes to the array external to this class will result in
     *            a change in the dimensions of this lattice (this is useful in
     *            resizing the lattice to the dimensions of a phase,
     *            by giving a reference to phase.boundary().dimension() here).
     *            The dimensions array reference is final and cannot be changed .
     *            The lattice dimension D is given by dimensions.length.
     * @param siteFactory
     *            makes the sites of the lattice
     */
    public CellLattice(Vector dimensions, SiteFactory siteFactory) {
        super(dimensions.D(), siteFactory);
        cellSize = new double[D()];
        idx = new int[D()];
        this.dimensions = dimensions;
    }
    
    /**
     * Returns the cell in which the given point lies. An ArrayIndexOutOfBounds
     * exception is thrown if the point lies outside the bounds of the lattice
     * (less than zero or greater than dimensions). The dimension of the vector
     * is assumed to be consistent with the dimension D of the lattice (i.e.,
     * r.D() == this.D()) but this is not checked.
     */
    public Object site(Vector r) {
        assignIndex(r, idx);
        return site(idx);

    }
    
    public void assignIndex(Vector r, int[] idx) {
        for(int i=0; i<idx.length; i++) {
            idx[i] = (int)(size[i]*r.x(i)/dimensions.x(i));
        }        
    }
    
    /**
     * Returns the array that specifies the spatial dimensions of the lattice of cells.
     */
    public Vector getDimensions() {
        return dimensions;
    }
    
    public void setDimensions(Vector d) {
        dimensions.E(d);
    }

    /**
     * Returns the spatial dimensions of all cells in the lattice, such that a cell in the
     * lattice will have a size in each direction as given by the corresponding value
     * in the returned array.  Value is computed on-the-fly using current value of dimensions
     * and size.
     */
    public double[] getCellSize() {
        for(int i=0; i<cellSize.length; i++) {
            cellSize[i] = dimensions.x(i)/size[i];
        }
        return cellSize;
    }
    
    private final double[] cellSize;
    private final Vector dimensions;
    private final int[] idx;//a work array

    /**
     * Extends the SimpleLattice neighbor iterator to provide methods that
     * specify the neighbor range in terms of a distance, rather than an index range. 
     */
    public static class NeighborIterator extends RectangularLatticeNbrIteratorSquare {
        public NeighborIterator(int D) {
            super(D);
            idx = new int[D];
            previousDimensions = Space.makeVector(D);
            neighborDistance = new double[D];
            setRange(Default.ATOM_SIZE);
        }
        
        /**
         * Defines neighbors to be all cells that are within the given distance
         * from a central cell.  Specifically, the integer range in direction
         * j is 1 + (int)(neighborDistance/cellSize[j]). 
         */
        public void setRange(double neighborDistance) {
            for(int i=0; i<D; i++) {
                this.neighborDistance[i] = neighborDistance;
            }
            setRange(this.neighborDistance);
        }

        /**
         * Defines neighbors to be all cells that are within the given distance
         * from a central cell, possible with a different distance specified in 
         * each direction.  Specifically, the integer range in direction
         * j is 1 + (int)(size[j]*neighborDistance/dimensions[j]). 
         */
        public void setRange(double[] neighborDistance) {
            if(neighborDistance.length != D) throw new IllegalArgumentException("Dimension of neighborDistance array inconsistent with dimension of lattice");
            if(lattice == null) return;
            for(int i=0; i<D; i++) {
                idx[i] = 1+(int)(lattice.getSize()[i]*neighborDistance[i]/((CellLattice)lattice).dimensions.x(i));
                this.neighborDistance[i] = neighborDistance[i];
            }
            super.setRange(idx);
        }
        
        /**
         * Checks whether lattice dimensions have changed, and updates integer
         * neighbor range if they have. 
         */
        public void checkDimensions() {
            if(lattice == null) return;
            Vector currentDimensions = ((CellLattice)lattice).getDimensions();
            if(!previousDimensions.equals(currentDimensions)) {
                setRange(neighborDistance);
                previousDimensions.E(currentDimensions);
            }
        }
        
        private final int[] idx;//a work array
        private final Vector previousDimensions;
        private final double[] neighborDistance;
    }//end of NeighborIterator
}
