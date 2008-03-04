/*
 * History
 * Created on Oct 27, 2004 by kofke
 */
package etomica.action;

import etomica.api.IController;

/**
 * Convenience class used to define a ControllerAction. Implements all methods
 * of ControllerAction interface, except for actionPerformed.
 */
public abstract class ControllerActionAdapter implements ControllerAction, java.io.Serializable {

	/**
	 * @return Returns the controller on which this action will be performed.
	 */
	public IController getController() {
		return controller;
	}

	/**
	 * @param controller
	 *            The controller on which this action will be performed.
	 */
	public void setController(IController controller) {
		this.controller = controller;
	}

	protected IController controller;
}
