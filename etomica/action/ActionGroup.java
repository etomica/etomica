package etomica.action;

import etomica.Action;


/**
 * A set of Action instances grouped and performed in series 
 * as if a single action.  Actions are defined at construction.
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Feb 2, 2005 by kofke
 */
public class ActionGroup implements Action {

    /**
     * Defines group via the given array of actions.  Copy
     * of array is made and used internally.  Assigns an uninformative
     * default label to group.
     */
    public ActionGroup(Action[] actions) {
        this("Action group", actions);
    }
    
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

    private final Action[] actions;
    private String label;
}
