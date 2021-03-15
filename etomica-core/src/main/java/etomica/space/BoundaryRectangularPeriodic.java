/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.IndexIteratorSizable;

/**
 * Rectangular boundary that is periodic in every dimension.
 */
public class BoundaryRectangularPeriodic extends BoundaryRectangular {

    protected final Vector tempImage;

    /**
     * Constructs cubic boundary with the default box-size given by the Simulation.
     */
    public BoundaryRectangularPeriodic(Space _space) {
        this(_space, 10.0);
    }

    /**
     * Constructs cubic boundary for the given Space, with each edge of length boxSize.
     */
    public BoundaryRectangularPeriodic(Space _space, double boxSize) {
        super(_space, boxSize);
        tempImage = space.makeVector();
        // call updateDimensions again so dimensionsHalf is updated
        updateDimensions();
    }

    /**
     * Constructs rectangular boundary for the given Space, with each edge of length boxSize.
     */
    public BoundaryRectangularPeriodic(Space _space, double[] boxSize) {
        super(_space, boxSize);
        tempImage = space.makeVector();
        // call updateDimensions again so dimensionsHalf is updated
        updateDimensions();
    }

    @Override
    public void updateDimensions() {
        super.updateDimensions();
    }

    @Override
    public IndexIteratorSizable getIndexIterator() {
        return new IndexIteratorRectangular(space.D());
    }

    @Override
    public Vector centralImage(Vector r) {
        tempImage.E(r);
        nearestImage(tempImage);
        tempImage.ME(r);
        return tempImage;
    }

    @Override
    public void nearestImage(Vector dr) {
        dr.nearestImage(dimensions);
    }

    @Override
    public boolean getPeriodicity(int d) {
        return true;
    }
}
