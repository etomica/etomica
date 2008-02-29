package etomica.action;

import etomica.api.IBox;
import etomica.box.Box;

/**
 * Convenience class used to define a BoxAction. Implements all methods
 * of BoxAction interface, except for actionPerformed.
 */

public abstract class BoxActionAdapter implements BoxAction, java.io.Serializable {

	/**
	 * @return Returns the box on which this action will be performed.
	 */
	public IBox getBox() {
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
