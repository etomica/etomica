package etomica;
import etomica.units.Dimension;

/**
 * A Modulator object permits changes to be made to a property of another object.
 * The property or object may be unknown to the object that is using the modulator.
 * For example, a modulator may be constructed to adjust a certain property, and then
 * the modulator can be passed to a generic Device, thereby connecting the device to the property.
 * @author David Kofke
 */

public abstract class ModulatorAbstract implements java.io.Serializable, DatumSource {
    
    /**
     * A string describing the property manipulated by the modulator
     */
    protected String label = new String("");

	/**
	 * Returns the physical dimensions (e.g., mass, length, pressure, etc.) of the quantity being modulated
	 */
    public abstract Dimension getDimension();
    
    public abstract void setValue(double d);
    public abstract double getValue();
    
    /**
     * Same as getValue, to implement DatumSource interface.  Argument is ignored.
     */
    public double value(DataSource.ValueType dummy) {return getValue();}
    
	/**
	 * Accessor method for the modulator's label
	 */
	public String getLabel() {return label;}
	/**
	 * Accessor method for the modulator's label
	 */
	public void setLabel(String s) {label = s;}
}