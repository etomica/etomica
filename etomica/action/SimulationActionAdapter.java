package etomica.action;

import etomica.simulation.Simulation;

/**
 * Convenience class used to define a SimulationAction. Implements all methods
 * of SimulationAction interface, except for actionPerformed.
 */
public abstract class SimulationActionAdapter implements SimulationAction, java.io.Serializable {

	/**
	 * Constructs the Action with the given descriptive label.
	 */
	public SimulationActionAdapter(String label) {
		setLabel(label);
	}

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


	/**
	 * Returns a descriptive label for this action.
	 */
	public String getLabel() {
		return label;
	}

	/**
	 * Sets a descriptive label for this action. This might be referenced, for
	 * example, by a button invoking this action in a graphical interface.
	 */
	public void setLabel(String label) {
		this.label = label;
	}

	protected Simulation simulation;

	private String label;
}
