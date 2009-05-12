package etomica.api;

public interface IEventManager {
	
	public void fireEvent(IEvent e);

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