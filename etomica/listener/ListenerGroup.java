package etomica.listener;

import etomica.api.IListener;

public interface ListenerGroup extends IListener {

    /**
     * Removes the given oldAction from the group.  oldAction must currently be
     * contained by this group.
     */
    public boolean removeListener(IListener oldAction);

    /**
     * Adds the given newAction to this group.  This group should not already
     * contain newAction.
     */
    public void addListener(IListener newAction);
    
    /**
     * Returns all actions from this group.
     */
    public IListener[] getAllListeners();
}
