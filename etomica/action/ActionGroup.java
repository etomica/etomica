package etomica.action;

import etomica.utility.Arrays;


/**
 * A set of Action instances grouped and performed in series 
 * as if a single action.  Actions may be defined at construction
 * and/or added afterward.  Actions are performed in the order in
 * which they are specifed in the constructor and subsequently added.
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Feb 2, 2005 by kofke
 */
public class ActionGroup implements Action, java.io.Serializable {

    /**
     * Constructs an action group that holds no actions.
     */
    public ActionGroup() {
        this(new Action[0]);
    }
    
    /**
     * Defines group via the given array of actions.  Copy
     * of array is made and used internally.  Assigns an uninformative
     * default label to group.
     */
    public ActionGroup(Action[] actions) {
        this("Action group", actions);
    }
    
    /**
     * Defines group via the given array of action and with the given label.
     * Copy of action array is made and used internally.
     */
    public ActionGroup(String label, Action[] actions) {
        this.actions = (Action[])actions.clone();
        setLabel(label);
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
    public void removeAction(Action oldAction) {
        actions = (Action[])Arrays.removeObject(actions, oldAction);
    }
    
    /**
     * Returns a string describing the group of actions.
     */
    public String getLabel() {
        return label;
    }

    /**
     * @param label The label to set.
     */
    public void setLabel(String label) {
        this.label = label;
    }

    private Action[] actions;
    private String label;
}
