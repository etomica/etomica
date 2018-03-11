/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.IndexIteratorSizable;

/**
 * Boundary that is not periodic in any direction.  Volume is specified,
 * but nothing about boundary enforces the dimensions on the atoms.  This effect must
 * be introduced with a potential or other construct that is made consistent
 * with the boundary.
 */
public class BoundaryRectangularNonperiodic extends BoundaryRectangular {


    /**
     * Make a boundary with unit volume.
     */
    public BoundaryRectangularNonperiodic(Space space) {
        super(space, 1.0);//boolean elements will all be false
        zero = space.makeVector();
    }

    public BoundaryRectangularNonperiodic(Space space, double size) {
        super(space, size);
        zero = space.makeVector();
    }

    /**
     * Returns a vector with all elements zero.
     */
    public Vector centralImage(Vector r) {
        zero.E(0.0);
        return zero;
    }

    /**
     * Does nothing.
     */
    public void nearestImage(Vector dr) {
    }

    /**
     * Returns a zero-length vector.
     */
    public double[][] imageOrigins(int nShells) {
        return origins;
    }

    public IndexIteratorSizable getIndexIterator() {
        return new IndexIteratorRectangular(0);
    }

    public boolean getPeriodicity(int d) {
        return false;
    }

    private final Vector zero;
    protected final double[][] origins= new double[0][0];//cannot be static because several boxs may be using at once
    private static final long serialVersionUID = 1L;
}
