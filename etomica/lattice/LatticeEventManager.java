package etomica.lattice;

import etomica.util.EventManager;

public class LatticeEventManager extends EventManager {

    public LatticeEventManager() {
        super();
        // TODO Auto-generated constructor stub
    }

    public void fireEvent(LatticeEvent event) {
        for(EventManager.Linker link=first; link!=null; link=link.next) {
            ((LatticeListener)link.listener).actionPerformed(event);
        }
    }

    protected Class getListenerClass() {return LatticeListener.class;}
}
