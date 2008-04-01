package etomica.api;

import etomica.simulation.SimulationEvent;

public interface ISimulationEventManager extends IEventManager {

	public abstract void fireEvent(SimulationEvent event);

}