/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.util.Arrays;

/**
 * A set of Action instances grouped and performed in series 
 * as if a single action.  Actions may be defined at construction
 * and/or added afterward.  Actions are performed in the order in
 * which they are specifed in the constructor and subsequently added.
 *
 * @author David Kofke
 *
 */
public class ActionGroupSeries implements IAction, java.io.Serializable, ActionGroup {

    /**
     * Constructs an action group that holds no actions.
     */
    public ActionGroupSeries() {
        this(new IAction[0]);
    }
    
    /**
     * Defines group via the given array of actions.  Copy
     * of array is made and used internally.
     */
    public ActionGroupSeries(IAction[] actions) {
        this.actions = actions.clone();
    }

    /**
     * Invokes the actionPerformed method of all actions
     * in the method, in the order given by the array at construction.
     */
    public void actionPerformed() {
        for(int i=0; i<actions.length; i++) {
            actions[i].actionPerformed();
        }
    }
    
    /**
     * Adds the given action to the group.  No check is made of whether
     * action is already in group; it is added regardless.  
     * @param newAction
     */
    public void addAction(IAction newAction) {
        actions = (IAction[])Arrays.addObject(actions, newAction);
    }
    
    /**
     * Removes the given action from the group.  No warning or
     * error is given if action is not in the group already.
     */
    public boolean removeAction(IAction oldAction) {
        int num = actions.length;
        actions = (IAction[])Arrays.removeObject(actions, oldAction);
        return actions.length != num; 
    }
    
    public IAction[] getAllActions() {
    	return actions.clone();
    }
    
    private static final long serialVersionUID = 1L;
    private IAction[] actions;
}
