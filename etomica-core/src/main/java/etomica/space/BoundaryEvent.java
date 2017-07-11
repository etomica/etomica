/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

/**
 * Event that informs a listener about some change to a boundary.
 */
public class BoundaryEvent {

    protected Boundary boundary;
    
    public BoundaryEvent(Boundary boundary) {
        this.boundary = boundary;
    }

    /**
     * @return the boundary which has changed.
     */
    public Boundary getBoundary() {
        return boundary;
    }
}
