/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;


/**
 * Interface for boundary objects that describe the size and periodic nature
 * of the box boundaries. Each Box has its own instance of this class. It
 * may be referenced by a coordinate pair when computing distances between
 * atoms, or by a cell iterator when generating neighbor lists. It is also used
 * by objects that enforce periodic images.
 * 
 * The boundary is responsible for firing inflate events when the boundary
 * dimensions change.
 */
public interface IBoundary {

    /**
     * Sets the box that holds the IBoundary.  If no box holds the boundary,
     * the box should be set to null.
     */
    public void setBox(IBox newBox);

    /**
     * Returns the boundary's IBox.  Might be null if the boundary is not
     * associated with a box.
     */
    public IBox getBox();

    /**
     * @return the volume enclosed by the boundary
     */
    public double volume();

    /**
     * Determines the translation vector needed to apply a periodic-image
     * transformation that moves the given point to an image point within the
     * boundary (if it lies outside, in a direction subject to periodic
     * imaging).
     * 
     * @param r
     *            vector position of untransformed point; r is not changed by
     *            this method
     * @return the displacement that must be applied to r to move it to its
     *         central-image location
     */
    public IVector centralImage(IVector r);

    /**
     * The nearest image is the pair of atom images that are closest when all
     * periodic-boundary images are considered.
     *
     * If the vector passed to this method is the displacement vector between
     * two points, the vector will be transformed such that it corresponds to
     * the vector between the nearest image of those two points.
     */
    public void nearestImage(IVectorMutable dr);

    /**
     * Returns the length of the sides of a rectangular box oriented in the lab
     * frame and in which the boundary is inscribed.  For a rectangular
     * boundary, this is simply the length of the boundary in each direction.
     * Each element of the returned vector gives in that coordinate direction
     * the maximum distance from one point on the boundary to another.
     * 
     * Manipulation of this copy will not cause any change to the boundary's
     * dimensions.
     */
    public IVector getBoxSize();

    /**
     * Scales the boundary dimensions such that the boundary's would be
     * inscribed within a rectangle of the of the given size.  For a
     * rectangular boundary, this simply sets the boundary length in each
     * dimension.  Specific interpretation of the given values for
     * non-rectangular shapes depends on the subclass.
     */
    public void setBoxSize(IVector v);

    /**
     * Returns the vector that defines the edge of this boundary for the given
     * dimension.  All vectors returned by this method can be considered to
     * originate from one corner.
     */
    public IVector getEdgeVector(int d);

    /**
     * Returns true if the boundary is periodic in the given direction (as
     * defined by the getEdgeVector method).
     */
   public boolean getPeriodicity(int d);
   
   /**
    * 
    * @return Returns the center point (origin) of the boundary
    */
   public IVector getCenter();
   
   public IBoundaryEventManager getEventManager();
}