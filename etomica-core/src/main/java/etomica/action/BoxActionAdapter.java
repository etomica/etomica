/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.box.Box;

/**
 * Convenience class used to define a BoxAction. Implements all methods
 * of BoxAction interface, except for actionPerformed.
 */

public abstract class BoxActionAdapter implements BoxAction, java.io.Serializable {

	/**
	 * @return Returns the box on which this action will be performed.
	 */
	public Box getBox() {
		return box;
	}

	/**
	 * @param box
	 *            The box on which this action will be performed.
	 */
	public void setBox(Box box) {
		this.box = box;
	}

	protected Box box;
}
