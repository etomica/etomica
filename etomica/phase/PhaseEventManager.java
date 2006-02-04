package etomica.phase;

import etomica.util.EventManager;

public class PhaseEventManager extends EventManager {

    public PhaseEventManager() {
        super();
        // TODO Auto-generated constructor stub
    }

    public void fireEvent(PhaseEvent event) {
        for(EventManager.Linker link=first; link!=null; link=link.next) {
            ((PhaseListener)link.listener).actionPerformed(event);
        }
    }

    protected Class getListenerClass() {return PhaseListener.class;}
}
