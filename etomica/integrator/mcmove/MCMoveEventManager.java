package etomica.integrator.mcmove;

import etomica.api.IEvent;
import etomica.api.IListener;
import etomica.util.EventManager;

public class MCMoveEventManager extends EventManager {

    public MCMoveEventManager() {
        super();
    }

    public void fireEvent(IEvent event) {
        for(EventManager.Linker link=first; link!=null; link=link.next) {
            ((IListener)link.listener).actionPerformed(event);
        }
    }

}
