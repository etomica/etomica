package etomica.integrator.mcmove;

import etomica.util.EventManager;
import etomica.util.IEvent;
import etomica.util.IListener;

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
