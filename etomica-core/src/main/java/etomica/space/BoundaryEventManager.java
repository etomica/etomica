/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

import java.util.ArrayList;
import java.util.List;

/**
 * Manager of listeners for boundary events. This object maintains a list of listeners
 * which receive the events.
 */
public class BoundaryEventManager {

    private final List<BoundaryEventListener> boundaryListeners = new ArrayList<>();

    /**
     * Adds the given listener to this event manager.
     *
     * @param newListener the listener to be added
     */
    public synchronized void addListener(BoundaryEventListener newListener) {
        if(newListener == null) throw new NullPointerException("Cannot add null as a listener to Box");
        if (boundaryListeners.contains(newListener)) {
            throw new RuntimeException(newListener+" is already an interval action");
        }
        boundaryListeners.add(newListener);
    }


    /**
     * Removes the given listener from this event manager.
     *
     * @param listener the listener to be removed
     */
    public synchronized void removeListener(BoundaryEventListener listener) {
        boundaryListeners.remove(listener);
    }

    public void inflate(Boundary boundary) {
        BoundaryEvent event = new BoundaryEvent(boundary);
        for (BoundaryEventListener boundaryListener : boundaryListeners) {
            boundaryListener.boundaryInflate(event);
        }
    }
}
