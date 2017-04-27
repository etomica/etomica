/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;

/**
 * Event that informs a listener about some change to a boundary.  The listener
 * can determine the specific type of even from the class of the event object
 * and can use event methods to determine what happened to the boundary.
 */
public interface IBoundaryEvent {
    
    /**
     * Returns the boundary which has changed.
     */
    public IBoundary getBoundary();
}
