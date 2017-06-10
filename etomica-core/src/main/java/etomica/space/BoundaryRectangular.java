/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

import etomica.lattice.IndexIteratorRectangular;
import etomica.math.geometry.*;

/**
 * Boundary that is in the shape of a rectangular parallelepiped.
 * Periodicity in each direction is specified by subclass.
 */
public abstract class BoundaryRectangular extends Boundary {

    private static final long serialVersionUID = 1L;
    protected final Vector dimensions;
    protected final float[][] shift0 = new float[0][0];
    protected final Vector[] edgeVectors;
    private final IndexIteratorRectangular indexIterator;

    /**
     * Constructs cubic boundary of the given periodicity, using the space and default box-size
     * given by the Simulation.
     */
    public BoundaryRectangular(Space _space) {
        this(_space, 10.0);
    }

    /**
     * Constructs cubic boundary of the given periodicity with each edge of length boxSize
     */
    public BoundaryRectangular(Space space, double boxSize) {
        this(space, makeArray(space.D(), boxSize));
    }

    /**
     * Constructs rectangular boundary of the given periodicity with edges given by the
     * values in the array boxSize.  Length of arrays must equal dimension of space.
     */
    public BoundaryRectangular(Space space, double[] boxSize) {
        super(space, makeShape(space));
        dimensions = space.makeVector();
        dimensions.E(boxSize);

        indexIterator = new IndexIteratorRectangular(space.D());
        edgeVectors = new Vector[space.D()];
        for (int i = 0; i < space.D(); i++) {
            edgeVectors[i] = space.makeVector();
        }
        updateDimensions();
    }

    private static final double[] makeArray(int n, double d) {
        double[] array = new double[n];
        for (int i = 0; i < n; i++) {
            array[i] = d;
        }
        return array;
    }

    //used by constructors
    private static Polytope makeShape(Space space) {
        switch (space.D()) {
            case 1:
                return new LineSegment(space);
            case 2:
                return new Rectangle(space);
            case 3:
                return new Cuboid(space);
            default:
                throw new IllegalArgumentException("BoundaryRectangular not appropriate to given space");
        }
    }

    /**
     * Returns a vector with elements equal to the lengths of the edges of
     * the boundary.  The returned Vector does not represent the values internally,
     * so manipulation of the vector has no effect on this BoundaryRectangular instance.
     */
    public Vector getBoxSize() {
        return dimensions;
    }

    /**
     * Sets the size and shape of the rectangular boundary.  Values are
     * copied, so manipulation of the given vector has no subsequent effect
     * on this Boundary instance.
     */
    public void setBoxSize(Vector v) {
        dimensions.E(v);
        updateDimensions();

        eventManager.inflate(this);
    }

    protected void updateDimensions() {
        ((Rectangular) shape).setEdgeLengths(dimensions);
        for (int i = 0; i < space.D(); i++) {
            edgeVectors[i].setX(i, dimensions.getX(i));
        }
    }

    /**
     * Returns the "volume" of the rectangular region defined by this Boundary.
     * For a 2D and 1D spaces, this volume is actually an area and length, respectively.
     */
    public double volume() {
        return shape.getVolume();
    }

    public Vector getEdgeVector(int d) {
        return edgeVectors[d];
    }

    /**
     * Returns a set of image origins for a set of periodic image shells.
     * The returned array is of dimension [(2*nShells+1)^D][D], where D
     * is the dimension of the space.
     */
    public double[][] imageOrigins(int nShells) {
        Vector workVector = space.makeVector();
        int shellFormula = (2 * nShells) + 1;
        int nImages = space.powerD(shellFormula) - 1;
        double[][] origins = new double[nImages][space.D()];
        indexIterator.setSize(shellFormula);
        indexIterator.reset();
        int k = 0;
        while (indexIterator.hasNext()) {
            int[] index = indexIterator.next();
            for (int i = 0; i < space.D(); i++) {
                workVector.setX(i, index[i]);
            }
            workVector.PE(-(double) nShells);
            if (workVector.isZero()) continue;
            workVector.TE(dimensions);
            workVector.assignTo(origins[k++]);
        }
        return origins;
    }
}
