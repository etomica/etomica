package etomica.units;

import java.io.ObjectStreamException;

/**
 * Unit class for the number, or quantity, of something. This is one of the
 * default simulation units.
 */
public final class Count extends SimpleUnit {

    /**
     * Singleton instance of this unit
     */
    public static final Count UNIT = new Count();

    private Count() {
        super(Quantity.DIMENSION, 1.0, "Count", "", Prefix.NOT_ALLOWED);
    }

    /**
     * Required to guarantee singleton when deserializing.
     * 
     * @return the singleton UNIT
     */
    private Object readResolve() throws ObjectStreamException {
        return UNIT;
    }
    
    private static final long serialVersionUID = 1;

}