package etomica.api;

public interface IEventManager {

	/**
	 * Adds a listener.  Synchronized to avoid conflict with removeListener.
	 */
	public abstract void addListener(Object listener);

	public abstract void addListener(Object listener, boolean doSerialize);

	/**
	 * Removes a listener.  Synchronized to avoid conflict with addListener.
	 */
	public abstract void removeListener(Object listener);

}