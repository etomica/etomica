package etomica.action.activity;

import etomica.util.EventManager;

public class ControllerEventManager extends EventManager {

    public ControllerEventManager() {
        super();
        // TODO Auto-generated constructor stub
    }

    public void fireEvent(ControllerEvent event) {
        for(EventManager.Linker link=first; link!=null; link=link.next) {
            ((ControllerListener)link.listener).actionPerformed(event);
        }
    }

    protected Class getListenerClass() {return ControllerListener.class;}
}
