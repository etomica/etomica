package etomica.simulation;

import etomica.util.EventManager;
import etomica.util.EventManager.Linker;

public class SimulationEventManager extends EventManager {

    public SimulationEventManager() {
        super();
        // TODO Auto-generated constructor stub
    }

    public void fireEvent(SimulationEvent event) {
        for(EventManager.Linker link=first; link!=null; link=link.next) {
            ((SimulationListener)link.listener).actionPerformed(event);
        }
    }

    protected Class getListenerClass() {return SimulationListener.class;}
}
