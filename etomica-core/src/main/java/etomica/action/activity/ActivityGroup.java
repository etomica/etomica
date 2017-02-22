/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action.activity;

import etomica.action.ActionGroup;
import etomica.action.IAction;

public interface ActivityGroup extends ActionGroup {

    /**
     * Returns all actions from this group that have been completed.
     */
    public IAction[] getCompletedActions();
    
    /**
     * Returns all actions from this group that are currently being performed.
     */
    public IAction[] getCurrentActions();
    
    /**
     * Returns all actions from this group that have not yet started.
     */
    public IAction[] getPendingActions();
}
