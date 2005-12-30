package etomica.units;

import java.io.ObjectStreamException;

/**
 * Simulation unit of time.
 */
public final class Picosecond extends SimpleUnit {

    /**
     * Singleton instance of this unit.
     */
    public static final Picosecond UNIT = new Picosecond();

    private Picosecond() {
        super(Time.DIMENSION, 1.0, "picoseconds", "ps", Prefix.NOT_ALLOWED);
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