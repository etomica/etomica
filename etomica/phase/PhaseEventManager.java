package etomica.phase;

import etomica.util.EventManager;

public class PhaseEventManager extends EventManager {

    public PhaseEventManager() {
        super();
    }

    public void fireEvent(PhaseEvent event) {
        for(EventManager.Linker link=first; link!=null; link=link.next) {
            ((PhaseListener)link.listener).actionPerformed(event);
        }
    }

    private static final long serialVersionUID = 1L;
    protected Class getListenerClass() {return PhaseListener.class;}
}
