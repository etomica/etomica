/*
 * History
 * Created on Oct 27, 2004 by kofke
 */
package etomica.action;

import etomica.api.IIntegrator;

/**
 * Convenience class used to define a IntegratorAction. Implements all methods
 * of IntegratorAction interface, except for actionPerformed.
 */
public abstract class IntegratorActionAdapter implements IntegratorAction, java.io.Serializable {

    public IntegratorActionAdapter() {}

    /**
	 * Constructs the Action with the given descriptive label.
	 */
	public IntegratorActionAdapter(IIntegrator integrator) {
	    setIntegrator(integrator);
	}

	/**
	 * @return Returns the integrator on which this action will be performed.
	 */
	public IIntegrator getIntegrator() {
		return integrator;
	}

	/**
	 * @param integrator
	 *            The integrator on which this action will be performed.
	 */
	public void setIntegrator(IIntegrator integrator) {
		this.integrator = integrator;
	}

	protected IIntegrator integrator;
}
