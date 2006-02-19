package etomica.integrator.mcmove;

import etomica.util.EventManager;

public class MCMoveEventManager extends EventManager {

    public MCMoveEventManager() {
        super();
    }

    public void fireEvent(MCMoveEvent event) {
        for(EventManager.Linker link=first; link!=null; link=link.next) {
            ((MCMoveListener)link.listener).actionPerformed(event);
        }
    }

    protected Class getListenerClass() {return MCMoveListener.class;}
}
