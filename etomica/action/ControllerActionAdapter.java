/*
 * History
 * Created on Oct 27, 2004 by kofke
 */
package etomica.action;

import etomica.Controller;

/**
 * Convenience class used to define a ControllerAction. Implements all methods
 * of ControllerAction interface, except for actionPerformed.
 */
public abstract class ControllerActionAdapter implements ControllerAction {

	/**
	 * Constructs the Action with the given descriptive label.
	 */
	public ControllerActionAdapter(String label) {
		setLabel(label);
	}

	/**
	 * @return Returns the controller on which this action will be performed.
	 */
	public Controller getController() {
		return controller;
	}

	/**
	 * @param controller
	 *            The controller on which this action will be performed.
	 */
	public void setController(Controller controller) {
		this.controller = controller;
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

	protected Controller controller;

}