/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

import etomica.api.IRandom;
import etomica.api.IVector;

/**
 * Interface for a class that specifies an orientation in space.
 * 
 * @author Andrew Schultz
 */
public interface IOrientation {

    /**
     * Copies the given orientation to this one.
     */
    void E(IOrientation o);

    /**
     * Returns a unit vector pointing in the orientation's direction.  This
     * vector should not be modified.
     */
    IVector getDirection();

    /**
     * Sets this orientation to point in the given direction.
     * @throws Exception if vector has 0 length
     */
    void setDirection(IVector newDirection);

    /**
     * Perform a rotation by a random amount in the solid angle theta on the 
     * present orientation.
     */
    void randomRotation(IRandom random, double theta);
}
