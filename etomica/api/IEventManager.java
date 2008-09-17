package etomica.api;

public interface IEventManager {

    public void fireEvent(IEvent event);

	public void addListener(IListener listener);

	public void addListener(IListener listener, boolean doSerialize);

	public void removeListener(IListener listener);

}