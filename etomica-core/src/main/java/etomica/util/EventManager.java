/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.util;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CopyOnWriteArrayList;

/**
 * Class to take care of listener lists and event firing for simulation elements.
 * A class can make an instance of this manager as a field, and delegate any
 * listener management functions to it.
 */
public class EventManager<E extends IEvent> {
    private final List<IListener<E>> listeners;

    public EventManager() {
        this.listeners = new ArrayList<>();
    }

    public void fireEvent(E e) {
        for(IListener<E> listener : listeners) {
            listener.actionPerformed(e);
        }
    }

    /**
     * Adds the given listener to this event manager.
     */
    public void addListener(IListener<E> listener) {
        listeners.add(listener);
    }

    /**
     * Removes the given listener from this event manager.
     */
    public void removeListener(IListener<E> listener) {
        this.listeners.remove(listener);

    }

}
