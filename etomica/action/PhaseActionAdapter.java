package etomica.action;

import etomica.Phase;

/**
 * Convenience class used to define a PhaseAction. Implements all methods
 * of PhaseAction interface, except for actionPerformed.
 */

public abstract class PhaseActionAdapter implements PhaseAction, java.io.Serializable {

	/**
	 * Constructs the Action with the given descriptive label.
	 */
	public PhaseActionAdapter(String label) {
		setLabel(label);
	}

	/**
	 * @return Returns the phase on which this action will be performed.
	 */
	public Phase getPhase() {
		return phase;
	}

	/**
	 * @param phase
	 *            The phase on which this action will be performed.
	 */
	public void setPhase(Phase phase) {
		this.phase = phase;
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

	private String label;

	protected Phase phase;

}//end of PhaseActionAdapter
