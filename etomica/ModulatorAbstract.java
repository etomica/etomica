package simulate;
import simulate.units.Dimension;

/**
 * A Modulator object permits changes to be made to a property of another object.
 * The property or object may be unknown to the object that is using the modulator.
 * For example, a modulator may be constructed to adjust a certain property, and then
 * the modulator can be passed to a generic Device, thereby connecting the device to the property.
 * Direct subclasses of ModulatorAbstract handle scalar and vector quantities, respectively.
 */

public abstract class ModulatorAbstract implements java.io.Serializable {
    
    /**
     * A string describing the property manipulated by the modulator
     */
    protected String label;

	/**
	 * Returns the physical dimensions (e.g., mass, length, pressure, etc.) of the quantity being modulated
	 */
    public abstract Dimension getDimension();
    
	/**
	 * Accessor method for the modulator's label
	 */
	public String getLabel() {return label;}
	/**
	 * Accessor method for the modulator's label
	 */
	public void setLabel(String s) {label = s;}
}