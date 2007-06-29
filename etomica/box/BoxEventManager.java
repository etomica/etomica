package etomica.box;

import etomica.util.EventManager;

public class BoxEventManager extends EventManager {

    public BoxEventManager() {
        super();
    }

    public void fireEvent(BoxEvent event) {
        for(EventManager.Linker link=first; link!=null; link=link.next) {
            ((BoxListener)link.listener).actionPerformed(event);
        }
    }

    private static final long serialVersionUID = 1L;
    protected Class getListenerClass() {return BoxListener.class;}
}
