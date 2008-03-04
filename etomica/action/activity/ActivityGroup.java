package etomica.action.activity;

import etomica.action.ActionGroup;
import etomica.api.IAction;

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
