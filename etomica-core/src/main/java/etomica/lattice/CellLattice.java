/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice;

import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.Debug;

/**
 * A lattice of arbitrarily-sized rectangular cells, such that a point in space can be 
 * associated with one of the cells.
 */
public class CellLattice extends RectangularLattice {

    /**
     * @param dimVector
     *            the spatial dimensions of the lattice, such that the total
     *            space occupied by the lattice cells a rectangular box with
     *            each side of length dimensions.x(i)/size[i] for side i. The
     *            dimensions at any time are given by the array passed here, so
     *            any changes to the array external to this class will result in
     *            a change in the dimensions of this lattice (this is useful in
     *            resizing the lattice to the dimensions of a box,
     *            by giving a reference to box.boundary().dimension() here).
     *            The dimensions array reference is final and cannot be changed .
     *            The lattice dimension D is given by dimensions.length.
     * @param siteFactory
     *            makes the sites of the lattice
     */
    public CellLattice(Space space, Vector dimVector, SiteFactory siteFactory) {
        super(space.D(), siteFactory);
        cellSize = new double[D()];
        this.dimensions = space.makeVector();
        dimensions.E(dimVector);
    }
    
    /**
     * Returns the cell in which the given point lies. An ArrayIndexOutOfBounds
     * exception is thrown if the point lies outside the bounds of the lattice
     * (less than zero or greater than dimensions). The dimension of the vector
     * is assumed to be consistent with the dimension D of the lattice (i.e.,
     * r.D() == this.D()) but this is not checked.
     */
    public Object site(Vector r) {
        int idx1D = 0;
        for(int i=0; i<D; i++) {
            if (Debug.ON && 2 * Math.abs(r.getX(i)) > dimensions.getX(i)) {
                throw new RuntimeException("vector " + r + " is outside the boundary " + dimensions);
            }
            int j = ((int)(size[i]*(r.getX(i)/dimensions.getX(i)+0.5)));
            if (j == -1) j = 0;
            else if (j == size[i]) j = size[i]-1;
            idx1D += j * jumpCount[i];
        }
        return sites[idx1D];
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
            cellSize[i] = dimensions.getX(i)/size[i];
        }
        return cellSize;
    }
    
    private static final long serialVersionUID = 1L;
    private final double[] cellSize;
    private final Vector dimensions;

    /**
     * Extends the SimpleLattice neighbor iterator to provide methods that
     * specify the neighbor range in terms of a distance, rather than an index range. 
     */
    public static class NeighborIterator extends RectangularLatticeNbrIteratorSquare {

        public NeighborIterator(int D, double range) {
            super(D);
            idx = new int[D];
            previousDimensions = Space.makeVector(D);
            setNeighborDistance(range);
        }

        public void setLattice(FiniteLattice lattice) {
            super.setLattice(lattice);
            // we got a lattice.  Perhaps setNeighborDistance was called before when we didn't
            setNeighborDistance(neighborDistance);
        }

        /**
         * Defines neighbors to be all cells that are within the given distance
         * from a central cell.  Specifically, the integer range in direction
         * j is 1 + (int)(neighborDistance/cellSize[j]). 
         */
        public void setNeighborDistance(double newNeighborDistance) {
            neighborDistance = newNeighborDistance;
            if(lattice == null) return;
            for(int i=0; i<D; i++) {
                double boundaryLength = ((CellLattice)lattice).getDimensions().getX(i);
                if (boundaryLength == 0) {
                    return ;
                }
                idx[i] = (int)(lattice.getSize()[i]*neighborDistance/boundaryLength);
                if (idx[i] < (lattice.getSize()[i]*neighborDistance/boundaryLength)) {
                    // if we rounded down, increment by 1 so we round up.
                    idx[i]++;
                }
                if (idx[i]*2 + 1 > lattice.getSize()[i]) {
                    // this can happen if the box is big enough for the potential cutoff,
                    // but not big enough to hold a lattice of cells of size
                    // neighborDistance/cellRange
                    if (lattice.getSize()[i] % 2 == 0) {
                        throw new RuntimeException("I really don't mind cheating a bit, but I need an odd number of cells.");
                    }
                    idx[i] = (lattice.getSize()[i]-1)/2;
                }
            }
            if (idx[0] == 0) throw new RuntimeException("not very nice");
            super.setRange(idx);
        }
        
        public double getNeighborDistance() {
            return neighborDistance;
        }
        
        public void reset() {
            checkDimensions();
            super.reset();
        }

        /**
         * Checks whether lattice dimensions have changed, and updates integer
         * neighbor range if they have. 
         */
        public void checkDimensions() {
            if(lattice == null) return;
            Vector currentDimensions = ((CellLattice)lattice).getDimensions();
            if(!previousDimensions.equals(currentDimensions)) {
                setNeighborDistance(neighborDistance);
                previousDimensions.E(currentDimensions);
            }
        }
        
        private static final long serialVersionUID = 1L;
        private final int[] idx;//a work array
        private final Vector previousDimensions;
        private double neighborDistance;
    }//end of NeighborIterator
}
