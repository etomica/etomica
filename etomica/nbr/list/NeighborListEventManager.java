/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.list;

import java.util.LinkedList;

public class NeighborListEventManager {

    private final LinkedList<INeighborListListener> intervalListeners = new LinkedList<INeighborListListener>();
    private static final long serialVersionUID = 1L;

    
    public synchronized void addListener(INeighborListListener newListener) {
        if(newListener == null) throw new NullPointerException("Cannot add null as a listener to Neighbor List");
        if (intervalListeners.contains(newListener)) {
            throw new RuntimeException(newListener+" is already an interval action");
        }
        intervalListeners.add(newListener);
    }

    public synchronized void removeListener(INeighborListListener listener) {
        intervalListeners.remove(listener);
    }

    public void neighborsUpdated() {
        for(int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).neighborListNeighborsUpdated();
        }
    }
}
