package etomica.simulation;

import etomica.api.IEvent;
import etomica.api.IListener;
import etomica.util.EventManager;

public class SimulationEventManager extends EventManager {

    public SimulationEventManager() {
        super();
    }

    /* (non-Javadoc)
	 * @see etomica.simulation.ISimulationEventManager#fireEvent(etomica.simulation.SimulationEvent)
	 */
    public void fireEvent(IEvent event) {
        for(EventManager.Linker link=first; link!=null; link=link.next) {
            ((IListener)link.listener).actionPerformed(event);
        }
    }

}
