package etomica.api;

public interface IBoxEventManager {

    /**
     * Adds the given listener to this event manager.
     */
    public void addListener(IBoxListener listener);

    /**
     * Removes the given listener from this event manager.
     */
    public void removeListener(IBoxListener listener);
}
