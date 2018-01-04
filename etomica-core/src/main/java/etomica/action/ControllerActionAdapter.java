/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * History
 * Created on Oct 27, 2004 by kofke
 */
package etomica.action;

import etomica.action.activity.Controller;

/**
 * Convenience class used to define a ControllerAction. Implements all methods
 * of ControllerAction interface, except for actionPerformed.
 */
public abstract class ControllerActionAdapter implements ControllerAction, java.io.Serializable {

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

	protected Controller controller;
}
