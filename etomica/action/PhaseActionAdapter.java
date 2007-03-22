package etomica.action;

import etomica.phase.Phase;

/**
 * Convenience class used to define a PhaseAction. Implements all methods
 * of PhaseAction interface, except for actionPerformed.
 */

public abstract class PhaseActionAdapter implements PhaseAction, java.io.Serializable {

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

	protected Phase phase;
}
