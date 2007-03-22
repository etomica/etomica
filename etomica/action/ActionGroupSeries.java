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
public class ActionGroupSeries implements Action, java.io.Serializable, ActionGroup {

    /**
     * Constructs an action group that holds no actions.
     */
    public ActionGroupSeries() {
        this(new Action[0]);
    }
    
    /**
     * Defines group via the given array of actions.  Copy
     * of array is made and used internally.
     */
    public ActionGroupSeries(Action[] actions) {
        this.actions = (Action[])actions.clone();
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
    public void addAction(Action newAction) {
        actions = (Action[])Arrays.addObject(actions, newAction);
    }
    
    /**
     * Removes the given action from the group.  No warning or
     * error is given if action is not in the group already.
     */
    public boolean removeAction(Action oldAction) {
        int num = actions.length;
        actions = (Action[])Arrays.removeObject(actions, oldAction);
        return actions.length != num; 
    }
    
    public Action[] getAllActions() {
    	return (Action[])actions.clone();
    }
    
    private static final long serialVersionUID = 1L;
    private Action[] actions;
}
