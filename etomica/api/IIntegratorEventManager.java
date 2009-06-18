package etomica.api;

public interface IIntegratorEventManager {

    /**
     * Adds the given listener to this event manager.
     */
    public void addListener(IIntegratorListener listener);

    /**
     * Removes the given listener from this event manager.
     */
    public void removeListener(IIntegratorListener listener);
    
}
