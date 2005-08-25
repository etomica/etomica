/*
 * History
 * Created on Oct 27, 2004 by kofke
 */
package etomica.action;

import etomica.integrator.Integrator;

/**
 * Convenience class used to define a IntegratorAction. Implements all methods
 * of IntegratorAction interface, except for actionPerformed.
 */
public abstract class IntegratorActionAdapter implements IntegratorAction, java.io.Serializable {

	/**
	 * Constructs the Action with the given descriptive label.
	 */
	public IntegratorActionAdapter(Integrator integrator, String label) {
	    setIntegrator(integrator);
        setLabel(label);
	}

	/**
	 * @return Returns the integrator on which this action will be performed.
	 */
	public Integrator getIntegrator() {
		return integrator;
	}

	/**
	 * @param integrator
	 *            The integrator on which this action will be performed.
	 */
	public void setIntegrator(Integrator integrator) {
		this.integrator = integrator;
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

	protected Integrator integrator;

}
