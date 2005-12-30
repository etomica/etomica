package etomica.units;

import java.io.ObjectStreamException;

/**
 * Null unit used for dimensionless quantities
 */
public final class Null extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Null();
    
    /**
     * Unit for an otherwise unspecified dimensionless quantity.
     */
    public static final Unit UNIT = new SimpleUnit(DIMENSION, 1.0, "", "", Prefix.NOT_ALLOWED);

    /**
     * Private constructor for singleton instantiation.
     */
    private Null() {
        super("Dimensionless", 0, 0, 0);
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