package etomica.api;

public interface IEventManager {

    /**
     * Notifies listeners held by this event manager of the given event.
     */
    public void fireEvent(IEvent event);

    /**
     * Adds the given listener to this event manager.
     */
	public void addListener(IListener listener);

	public void addListener(IListener listener, boolean doSerialize);

	/**
	 * Removes the given listener from this event manager.
	 */
	public void removeListener(IListener listener);

}