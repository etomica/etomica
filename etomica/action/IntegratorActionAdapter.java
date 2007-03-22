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

    public IntegratorActionAdapter() {}

    /**
	 * Constructs the Action with the given descriptive label.
	 */
	public IntegratorActionAdapter(Integrator integrator) {
	    setIntegrator(integrator);
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

	protected Integrator integrator;
}
