package etomica.box;

import etomica.api.IEvent;
import etomica.api.IListener;
import etomica.util.EventManager;

public class BoxEventManager extends EventManager {

    public BoxEventManager() {
        super();
    }

    /* (non-Javadoc)
	 * @see etomica.box.IBoxEventManager#fireEvent(etomica.box.BoxEvent)
	 */
    public void fireEvent(IEvent event) {
        for(EventManager.Linker link=first; link!=null; link=link.next) {
            ((IListener)link.listener).actionPerformed(event);
        }
    }

    private static final long serialVersionUID = 1L;

}
