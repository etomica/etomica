/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space3d;

import etomica.space.Vector;

/**
 * Interface representing an orientation in 3D.  The orientation defined by two
 * unit vectors that are orthogonal to one another.  They define two axes of
 * the orientation, with the third given by their cross product.
 * 
 * @author Andrew Schultz
 */
public interface IOrientationFull3D extends IOrientation3D {

    /**
     * Returns a unit vector pointing in the orientation's secondary direction.
     * This vector should not be modified.
     */
    public Vector getSecondaryDirection();
    
    /**
     * Sets the orientation's primary and secondary direction to be the given
     * directions.  The two vectors should be orthogonal.
     */
    public void setDirections(Vector newPrimaryDirection, Vector newSecondaryDirection);
}
