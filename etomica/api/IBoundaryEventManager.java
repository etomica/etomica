package etomica.api;

public interface IBoundaryEventManager {

    /**
     * Adds the given listener to this event manager.
     */
    public void addListener(IBoundaryListener listener);

    public void addListener(IBoundaryListener listener, boolean doSerialize);

    /**
     * Removes the given listener from this event manager.
     */
    public void removeListener(IBoundaryListener listener);
}
