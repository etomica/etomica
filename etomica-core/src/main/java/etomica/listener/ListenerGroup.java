/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.listener;

import etomica.util.IListener;

public interface ListenerGroup extends IListener {

    /**
     * Removes the given oldAction from the group.  oldAction must currently be
     * contained by this group.
     */
    public boolean removeListener(IListener oldAction);

    /**
     * Adds the given newAction to this group.  This group should not already
     * contain newAction.
     */
    public void addListener(IListener newAction);
    
    /**
     * Returns all actions from this group.
     */
    public IListener[] getAllListeners();
}
