/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

import etomica.util.random.IRandom;

/**
 * Interface for a class that specifies an orientation in space.
 * 
 * @author Andrew Schultz
 */
public interface IOrientation {

    /**
     * Copies the given orientation to this one.
     *
     * @param o the new orientation
     */
    void E(IOrientation o);

    /**
     * @ return a unit vector pointing in the orientation's direction.  This
     * vector should not be modified.
     */
    Vector getDirection();

    /**
     * Sets this orientation to point in the given direction.
     *
     * @param newDirection the new direction for this orientation
     * @throws Exception if vector has 0 length
     */
    void setDirection(Vector newDirection);

    /**
     * Perform a rotation by a random amount in the solid angle theta on the 
     * present orientation.
     *
     * @param random random number generator used to select the angle
     * @param theta the maximum angle (in radians) for random rotation
     */
    void randomRotation(IRandom random, double theta);
}
