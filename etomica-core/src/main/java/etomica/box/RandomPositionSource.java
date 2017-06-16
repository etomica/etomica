/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.space.Vector;

/**
 * Interface for an object which returns random positions within a box's
 * boundary.
 *
 * @author Andrew Schultz
 */
public interface RandomPositionSource {

    /**
     * Notifies the RandomPositionSource of the box from which random positions
     * should be taken from.
     */
    void setBox(Box box);

    /**
     * Returns a random position with the previously set box's boundary.
     */
    Vector randomPosition();
}
