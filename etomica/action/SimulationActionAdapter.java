package etomica.action;

import etomica.simulation.Simulation;

/**
 * Convenience class used to define a SimulationAction. Implements all methods
 * of SimulationAction interface, except for actionPerformed.
 */
public abstract class SimulationActionAdapter implements SimulationAction, java.io.Serializable {

	/**
	 * @return Returns the simulation on which this action will be performed.
	 */
	public Simulation getSimulation() {
		return simulation;
	}

	/**
	 * @param simulation
	 *            The simulation on which this action will be performed.
	 */
	public void setSimulation(Simulation simulation) {
		this.simulation = simulation;
	}

	protected Simulation simulation;
}
