package etomica.api;

import etomica.box.BoxEvent;

public interface IBoxEventManager {

	public abstract void fireEvent(BoxEvent event);
	
	public void addListener(Object listener);
	
	public void addListener(Object listener, boolean doSerialize);
	
	public void removeListener(Object listener);

}