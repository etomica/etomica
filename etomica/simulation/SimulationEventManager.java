package etomica.simulation;

import etomica.api.ISimulationEventManager;
import etomica.util.EventManager;

public class SimulationEventManager extends EventManager implements ISimulationEventManager {

    public SimulationEventManager() {
        super();
    }

    /* (non-Javadoc)
	 * @see etomica.simulation.ISimulationEventManager#fireEvent(etomica.simulation.SimulationEvent)
	 */
    public void fireEvent(SimulationEvent event) {
        for(EventManager.Linker link=first; link!=null; link=link.next) {
            ((SimulationListener)link.listener).actionPerformed(event);
        }
    }

    protected Class getListenerClass() {return SimulationListener.class;}
}
