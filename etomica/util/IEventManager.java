/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.util;


public interface IEventManager {
	
	public void fireEvent(IEvent e);

    /**
     * Adds the given listener to this event manager.
     */
	public void addListener(IListener listener);

	/**
	 * Removes the given listener from this event manager.
	 */
	public void removeListener(IListener listener);
}