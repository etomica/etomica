package etomica.units;

import java.io.ObjectStreamException;

/**
 * Undefined dimension used for quantities with undefined or unknown dimensions.
 */
public final class Undefined extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Undefined();
    /**
     * Unit for a quantity having undefined or unknown units. 
     * Any conversion using this unit will cause output to be NaN (not a number).
     */
    public static final Unit UNIT = new SimpleUnit(DIMENSION, Double.NaN, "undefined", "", Prefix.NOT_ALLOWED);

    /**
     * Private constructor for singleton instantiation.
     */
    private Undefined() {
        super("Undefined",
                Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN);
    }
    
    /**
     * Returns true if the given object is the singleton instance of this class.
     */
    public boolean equals(Object obj) {
        return obj == this;
    }

    /**
     * Returns UNIT, regardless of the given unit system.
     */
    public Unit getUnit(UnitSystem unitSystem) {
        return UNIT;
    }
    
    /**
     * Required to guarantee singleton when deserializing.
     * 
     * @return the singleton DIMENSION
     */
    private Object readResolve() throws ObjectStreamException {
        return DIMENSION;
    }

    private static final long serialVersionUID = 1;
}