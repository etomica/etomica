/*
 * Created on Jan 29, 2004
 */
package etomica.chem;
import java.awt.Color;

/**
 * Abstract structure for a class defining one of the chemical elements.
 * Subclasses are (or will be) defined to correspond to each element in the
 * periodic table.
 */
public abstract class Element {
	private Color color;
	
	public Element(Color color) {
		this.color = color;
	}
	/**
	 * Returns the color used when displaying the atom graphically.
	 * @return java.awt.Color
	 */
	public Color getColor() {
		return color;
	}

	/**
	 * Sets the color.
	 * @param color The color to set
	 */
	public void setColor(Color color) {
		this.color = color;
	}

}
