/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;

/**
 * Manager for box events.  This object maintains a list of listeners which
 * receive the events.
 */
public interface IBoxEventManager {

    /**
     * Adds the given listener to this event manager.
     *
     * @param listener the listener to be added
     */
    public void addListener(IBoxListener listener);

    /**
     * Removes the given listener from this event manager.
     *
     * @param listener the listener to be added
     */
    public void removeListener(IBoxListener listener);
}
