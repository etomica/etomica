/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;

import etomica.space.Boundary;

/**
 * Event that informs a listener about some change to a boundary.
 */
public interface IBoundaryEvent {
    
    /**
     * @return the boundary which has changed.
     */
    public Boundary getBoundary();
}
