/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;

import etomica.space.BoundaryEvent;

/**
 * Listener for boundary events.  The listener methods will be called when the
 * boundary changes.
 */
public interface IBoundaryListener {
    
    /**
     * Informs the listener that the boundary shape and/or size has changed.
     *
     * @param e event, which can be used to determine the boundary that changed.
     */
    public void boundaryInflate(BoundaryEvent e);
    
}
