/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.IndexIteratorSizable;

/**
 * Class for implementing slit periodic boundary conditions, in which
 * one dimension is not periodic.  Selection of non-periodic dimension may be changed
 * after construction.
 */
public class BoundaryRectangularSlit extends BoundaryRectangular {
    
    /**
     * Makes cubic volume with the x-dimension (index 0)
     * not periodic.  Length of each box edge is given by default boxSize in
     * given Simulation.
     */
    public BoundaryRectangularSlit(ISpace _space) {
        //consumer can set appropriate slit dim later
        this(0, _space);
    }
    
    /**
     * Makes cubic volume with the indicated dimension not periodic.
     * Length of each box edge is given by the default boxSize in the
     * given Simulation.
     * 
     * @param slitDim index indicating dimension that is not periodic (0 for x-dimension,
     * 1 for y-dimension, etc.).
     * @throws IllegalArgumentException if not (0 <= slitDim < space.D).
     */
    public BoundaryRectangularSlit(int slitDim, ISpace _space) {
        this(slitDim, 10.0, _space);
    }
    
    /**
     * Constructor for periodic boundary conditions with a slit 
     * in the given dimension.
     * @param space
     * @param slitDim slit dimension (in which PBC is not imposed).
     */
    public BoundaryRectangularSlit(int slitDim, double boxSize, ISpace _space) {
        super(_space,boxSize);
        sDim = slitDim;
        dimensionsHalf = space.makeVector();
        tempImage = space.makeVector();
        // call updateDimensions again so dimensionsHalf is updated
        updateDimensions();
    }

    /**
     * Sets the non-periodic dimension to that indicated by the given value 
     * (0 is x-dimension, 1 is y-dimension, etc.).
     * 
     * @throws IllegalArgumentException if not (0 <= slitDim < space.D).
     */
    public void setSlitDim(int slitDim) {
        sDim = slitDim;
    }
    
    /**
     * Returns index of non-periodic dimension.
     */
    public int getSlitDim() {
        return sDim;
    }
    
    protected void updateDimensions() {
        super.updateDimensions();
        // superclass constructor calls this before dimensionsHalf has been instantiated
        if (dimensionsHalf != null) {
            dimensionsHalf.Ea1Tv1(0.5,dimensions);
        }
    }

    public void nearestImage(IVectorMutable dr) {
        double x = dr.getX(sDim);
        dr.PE(dimensionsHalf);
        dr.mod(dimensions);
        dr.ME(dimensionsHalf);
        dr.setX(sDim,x);
    }
    
    public IVector centralImage(IVector r) {
        tempImage.E(r);
        nearestImage(tempImage);
        tempImage.ME(r);
        return tempImage;
    }

    public IndexIteratorSizable getIndexIterator() {
        return new IndexIteratorRectangular(space.D()-1);
    }

    public boolean getPeriodicity(int d) {
        return d != sDim;
    }

    private int sDim;
    private static final long serialVersionUID = 1L;
    protected final IVectorMutable dimensionsHalf;
    protected final IVectorMutable tempImage;
}
