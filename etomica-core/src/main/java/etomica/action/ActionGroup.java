/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;


public interface ActionGroup extends IAction {

    /**
     * Removes the given oldAction from the group.  oldAction must currently be
     * contained by this group.
     */
    public boolean removeAction(IAction oldAction);

    /**
     * Adds the given newAction to this group.  This group should not already
     * contain newAction.
     */
    public void addAction(IAction newAction);
    
    /**
     * Returns all actions from this group.
     */
    public IAction[] getAllActions();

}