/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.IndexIteratorSizable;

/**
 * Class for implementing pore periodic boundary conditions, in which
 * one dimension is periodic.  Selection of periodic dimension may be changed
 * after construction.
 */
public class BoundaryRectangularPore extends BoundaryRectangular {
    
    /**
     * Makes cubic volume with the x-dimension (index 0)
     * not periodic.  Length of each box edge is given by default boxSize in
     * given Simulation.
     */
    public BoundaryRectangularPore(Space space) {
        //consumer can set appropriate slit dim later
        this(space,0);
    }
    
    /**
     * Makes cubic volume with the indicated dimension not periodic.
     * Length of each box edge is given by the default boxSize in the
     * given Simulation.
     * 
     * @param poreDim index indicating dimension that is not periodic (0 for x-dimension,
     * 1 for y-dimension, etc.).
     * @throws IllegalArgumentException if not (0 <= slitDim < space.D).
     */
    public BoundaryRectangularPore(Space space, int poreDim) {
        this(space, poreDim, 10.0);
    }
    
    /**
     * Constructor for periodic boundary conditions with a pore
     * in the given dimension.
     * @param space
     * @param poreDim pore dimension (in which PBC is not imposed).
     */
    public BoundaryRectangularPore(Space space, int poreDim, double boxSize) {
        super(space,boxSize);
        pDim = poreDim;
        dimensionsHalf = space.makeVector();
        tempImage = space.makeVector();
        // call updateDimensions again so dimensionsHalf is updated
        updateDimensions();
    }

    /**
     * Sets the periodic dimension to that indicated by the given value 
     * (0 is x-dimension, 1 is y-dimension, etc.).
     * 
     * @throws IllegalArgumentException if not (0 <= slitDim < space.D).
     */
    public void setPoreDim(int poreDim) {
        pDim = poreDim;
    }
    
    /**
     * Returns index of periodic dimension.
     */
    public int getPoreDim() {
        return pDim;
    }
    
    protected void updateDimensions() {
        super.updateDimensions();
        // superclass constructor calls this before dimensionsHalf has been instantiated
        if (dimensionsHalf != null) {
            dimensionsHalf.Ea1Tv1(0.5,dimensions);
        }
    }

    public void nearestImage(Vector dr) {
        double x = dr.getX(pDim);
        if (x < -dimensionsHalf.getX(pDim)) {
            x += dimensions.getX(pDim);
            dr.setX(pDim, x);
        }
        else if (x > dimensionsHalf.getX(pDim)) {
            x -= dimensions.getX(pDim);
            dr.setX(pDim, x);
        }
    }
    
    public Vector centralImage(Vector r) {
        tempImage.E(r);
        nearestImage(tempImage);
        tempImage.ME(r);
        return tempImage;
    }


    public IndexIteratorSizable getIndexIterator() {
        return new IndexIteratorRectangular(1);
    }

    public boolean getPeriodicity(int d) {
        return d == pDim;
    }

    private int pDim;
    private static final long serialVersionUID = 1L;
    protected final Vector dimensionsHalf;
    protected final Vector tempImage;
}
