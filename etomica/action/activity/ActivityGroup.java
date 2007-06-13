package etomica.action.activity;

import etomica.action.Action;
import etomica.action.ActionGroup;

public interface ActivityGroup extends ActionGroup {

    /**
     * Returns all actions from this group that have been completed.
     */
    public Action[] getCompletedActions();
    
    /**
     * Returns all actions from this group that are currently being performed.
     */
    public Action[] getCurrentActions();
    
    /**
     * Returns all actions from this group that have not yet started.
     */
    public Action[] getPendingActions();
}
