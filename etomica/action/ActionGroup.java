package etomica.action;

public interface ActionGroup extends Action {

    /**
     * Removes the given oldAction from the group.  oldAction must currently be
     * contained by this group.
     */
    public boolean removeAction(Action oldAction);

    /**
     * Adds the given newAction to this group.  This group should not already
     * contain newAction.
     */
    public void addAction(Action newAction);
    
    /**
     * Returns all actions from this group.
     */
    public Action[] getAllActions();

}