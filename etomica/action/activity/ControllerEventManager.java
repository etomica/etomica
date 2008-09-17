package etomica.action.activity;

import etomica.api.IEvent;
import etomica.api.IListener;
import etomica.util.EventManager;

public class ControllerEventManager extends EventManager {

    public ControllerEventManager() {
        super();
    }

    /* (non-Javadoc)
	 * @see etomica.action.activity.IControllerEventManager#fireEvent(etomica.action.activity.ControllerEvent)
	 */
    public void fireEvent(IEvent event) {
        for(EventManager.Linker link=first; link!=null; link=link.next) {
            ((IListener)link.listener).actionPerformed(event);
        }
    }

//    protected Class getListenerClass() {return ControllerListener.class;}
}
