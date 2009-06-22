package etomica.action.activity;

import etomica.util.EventManager;
import etomica.util.IEvent;
import etomica.util.IListener;

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

}
