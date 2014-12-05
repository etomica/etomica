/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space2d;

import etomica.space.IOrientation;

/**
 * Interface for a class that specifies an orientation in a 2D space.
 */
public interface IOrientation2D extends IOrientation {
    
    /**
     * Rotate the orientation by the amount dt (radians) in the
     * counter-clockwise direction.
     */
    public void rotateBy(double dt);
}
