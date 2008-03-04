package etomica.action.activity;

import etomica.api.IControllerEventManager;
import etomica.util.EventManager;

public class ControllerEventManager extends EventManager implements IControllerEventManager {

    public ControllerEventManager() {
        super();
    }

    /* (non-Javadoc)
	 * @see etomica.action.activity.IControllerEventManager#fireEvent(etomica.action.activity.ControllerEvent)
	 */
    public void fireEvent(ControllerEvent event) {
        for(EventManager.Linker link=first; link!=null; link=link.next) {
            ((ControllerListener)link.listener).actionPerformed(event);
        }
    }

    protected Class getListenerClass() {return ControllerListener.class;}
}
