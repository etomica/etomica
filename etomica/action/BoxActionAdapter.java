package etomica.action;

import etomica.api.IBox;

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
	public void setBox(IBox box) {
		this.box = box;
	}

	protected IBox box;
}
