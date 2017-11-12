/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

/**
 * An action formed from a set of other actions.
 */
public interface ActionGroup extends IAction {

    /**
     * Removes a specified action from the group. If action appears multiple times
     * in group, first instance in execution order is removed.
     *
     * @param oldAction the action to be removed
     * @return false if the given action is not in the group; true otherwise
     */
    public boolean removeAction(IAction oldAction);

    /**
     * Adds the given action to this group, placing it at the end of the list.  If action is
     * already in list, it will be performed repeatedly.
     *
     * @param newAction the action to be added
     */
    public void addAction(IAction newAction);
    
    /**
     * Returns all actions from this group. List of actions held by group is unchanged.
     */
    public IAction[] getAllActions();

}