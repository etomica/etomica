package etomica.api;

import etomica.action.activity.ControllerEvent;

public interface IControllerEventManager {

	public abstract void fireEvent(ControllerEvent event);

    public void addListener(Object listener);
    
    public void addListener(Object listener, boolean doSerialize);
}