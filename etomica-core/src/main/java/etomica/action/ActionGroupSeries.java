/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import java.util.ArrayList;
import java.util.List;

/**
 * A set of Action instances grouped and performed in series
 * as if a single action.  Actions may be defined at construction
 * and/or added afterward.  Actions are performed in the order in
 * which they are specified in the constructor and subsequently added.
 *
 * @author David Kofke
 */
public class ActionGroupSeries implements ActionGroup {

    private final List<IAction> actions = new ArrayList<>();

    /**
     * Defines group via the given set of actions.
     */
    public ActionGroupSeries(IAction... actionArray) {
        for(IAction action : actionArray) {
            addAction(action);
        }
    }

    /**
     * Invokes all actions in the group, in the order that they were added at construction and
     * in previous calls to addAction.
     */
    public void actionPerformed() {
        for(IAction action : actions) action.actionPerformed();
    }

    public void addAction(IAction newAction) {
        actions.add(newAction);
    }

    public boolean removeAction(IAction oldAction) {
        return actions.remove(oldAction);
    }

    public IAction[] getAllActions() {
        return actions.toArray(new IAction[] {});
    }
}
