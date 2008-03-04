package etomica.box;

import etomica.api.IBoxEventManager;
import etomica.util.EventManager;

public class BoxEventManager extends EventManager implements IBoxEventManager {

    public BoxEventManager() {
        super();
    }

    /* (non-Javadoc)
	 * @see etomica.box.IBoxEventManager#fireEvent(etomica.box.BoxEvent)
	 */
    public void fireEvent(BoxEvent event) {
        for(EventManager.Linker link=first; link!=null; link=link.next) {
            ((BoxListener)link.listener).actionPerformed(event);
        }
    }

    private static final long serialVersionUID = 1L;
    protected Class getListenerClass() {return BoxListener.class;}
}
