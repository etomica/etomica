package etomica.api;

import etomica.simulation.SimulationEvent;

public interface ISimulationEventManager {

	public abstract void fireEvent(SimulationEvent event);

}