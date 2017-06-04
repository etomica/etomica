/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

import etomica.api.IBoundaryEvent;
import etomica.api.IBoundaryEventManager;
import etomica.box.Box;
import etomica.lattice.IndexIteratorSizable;
import etomica.math.geometry.Polytope;

/**
 * Parent class of boundary objects that describe the size and periodic nature
 * of the box boundaries. Each Box has its own instance of this class. It
 * may be referenced by a coordinate pair when computing distances between
 * atoms, or by a cell iterator when generating neighbor lists. It is also used
 * by objects that enforce periodic images.
 * <p>
 * The boundary is responsible for firing inflate events when the boundary
 * dimensions change.
 */
public abstract class Boundary {

    protected final Polytope shape;
    protected final Space space;
    private final Vector center;
    protected Box box;
    protected IBoundaryEvent inflateEvent;
    protected BoundaryEventManager eventManager;

    /**
     * Subclasses must invoke this constructor and provide a Space instance that
     * can be used to generate Vectors, and a Polytope that defines the shape
     * and volume of the boundary region. Both are final.
     */
    public Boundary(Space space, Polytope shape) {
        this.space = space;
        this.shape = shape;
        double zip[] = new double[space.D()];
        for (int i = 0; i < space.D(); i++) zip[i] = 0.0;
        center = space.makeVector(zip);
        eventManager = new BoundaryEventManager();
    }

    /**
     * @return the boundary's Box.  Might be null if the boundary is not
     * associated with a box.
     */
    public Box getBox() {
        return box;
    }

    /**
     * Sets the box that holds the IBoundary.  If no box holds the boundary,
     * the box should be set to null.
     *
     * @param newBox the box that holds the boundary
     */
    public void setBox(Box newBox) {
        box = newBox;
        if (box == null) {
            inflateEvent = null;
        } else {
            inflateEvent = new BoundaryEvent(this);
        }
    }

    public Polytope getShape() {
        return shape;
    }

    /**
     * @return the volume enclosed by the boundary
     */
    public double volume() {
        return shape.getVolume();
    }

    /**
     * @return the center point (origin) of the boundary
     */
    public Vector getCenter() {
        return center;
    }

    /**
     * @return the event manager, which fires notifications about changes to
     * this boundary to any added listener.
     */
    public IBoundaryEventManager getEventManager() {
        return eventManager;
    }

    public abstract IndexIteratorSizable getIndexIterator();

    /**
     * Set of vectors describing the displacements needed to translate the
     * central image to all of the periodic images. The first index specifies
     * each periodic image, while the second index indicates the xyz components
     * of the translation vector.
     *
     * @param nShells the number of shells of images to be computed
     */
    public abstract double[][] imageOrigins(int nShells);

    /**
     * Determines the translation vector needed to apply a periodic-image
     * transformation that moves the given point to an image point within the
     * boundary (if it lies outside, in a direction subject to periodic
     * imaging).
     *
     * @param r vector position of untransformed point; r is not changed by
     *          this method
     * @return the displacement that must be applied to r to move it to its
     * central-image location
     */
    public abstract Vector centralImage(Vector r);

    /**
     * The nearest image is the pair of atom images that are closest when all
     * periodic-boundary images are considered.
     * <p>
     * If the vector passed to this method is the displacement vector between
     * two points, the vector will be transformed such that it corresponds to
     * the vector between the nearest image of those two points.
     *
     * @param dr the vector to be transformed
     */
    public abstract void nearestImage(Vector dr);

    /**
     * Returns the length of the sides of a rectangular box oriented in the lab
     * frame and in which the boundary is inscribed.  For a rectangular
     * boundary, this is simply the length of the boundary in each direction.
     * Each element of the returned vector gives in that coordinate direction
     * the maximum distance from one point on the boundary to another.
     * <p>
     * Manipulation of this copy will not cause any change to the boundary's
     * dimensions.
     *
     * @return the box size
     */
    public abstract Vector getBoxSize();

    /**
     * Scales the boundary dimensions such that the boundary's would be
     * inscribed within a rectangle of the of the given size.  For a
     * rectangular boundary, this simply sets the boundary length in each
     * dimension.  Specific interpretation of the given values for
     * non-rectangular shapes depends on the subclass.
     *
     * @param v the box's new size
     */
    public abstract void setBoxSize(Vector v);

    /**
     * Returns the vector that defines the edge of this boundary for the given
     * dimension.  All vectors returned by this method can be considered to
     * originate from one corner.
     *
     * @param d the dimension of the desired edge vector
     * @return the edge vector
     */
    public abstract Vector getEdgeVector(int d);

    /**
     * Returns true if the boundary is periodic in the given direction (as
     * defined by the getEdgeVector method).
     *
     * @param d the dimension of the desired periodicity
     * @return the periodicity of dimension d
     */
    public abstract boolean getPeriodicity(int d);
}
