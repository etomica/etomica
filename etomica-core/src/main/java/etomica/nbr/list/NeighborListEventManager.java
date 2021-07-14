/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.list;

import java.util.ArrayList;
import java.util.List;

public class NeighborListEventManager {

    private final List<INeighborListener> listeners = new ArrayList<>();

    public void addListener(INeighborListener newListener) {
        if (newListener == null) {
            throw new NullPointerException("Cannot add null as a listener to Neighbor List");
        }
        if (listeners.contains(newListener)) {
            throw new RuntimeException(newListener + " is already an interval action");
        }
        listeners.add(newListener);
    }

    public void removeListener(INeighborListener listener) {
        listeners.remove(listener);
    }

    public void neighborsUpdated() {
        for (INeighborListener listener : listeners) {
            listener.neighborListNeighborsUpdated();
        }
    }
}
