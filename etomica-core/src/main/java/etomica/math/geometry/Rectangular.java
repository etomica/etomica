/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.geometry;

import etomica.space.Vector;


/**
 * Interface for a polytope with shape that can be specified by the length
 * of D sides.  For example, a rectangle in 2D, or a cuboid in 3D.
 *
 * @author David Kofke
 *
 */

public interface Rectangular {

    /**
     * Specifies the D edge lengths needed to describe the shape
     * of the polytope.  The given vector should be for a D-dimensional
     * space, and its values will be used to set the corresponding edge lengths.
     */
    public void setEdgeLengths(Vector v);
    
    /**
     * @return the D edge lengths as the elements of the vector.
     */
    public Vector getEdgeLengths();
}
