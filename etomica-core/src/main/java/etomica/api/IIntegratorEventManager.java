/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;

public interface IIntegratorEventManager {

    /**
     * Adds the given listener to this event manager.
     */
    public void addListener(IIntegratorListener listener);

    /**
     * Removes the given listener from this event manager.
     */
    public void removeListener(IIntegratorListener listener);

    /**
     * Returns true if the event manager is currently firing events.
     */
    public boolean firingEvent();
}
