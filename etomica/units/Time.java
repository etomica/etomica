package etomica.units;

import java.io.ObjectStreamException;

import etomica.units.systems.UnitSystem;

/**
 * Dimension for all units of time.
 */
public final class Time extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Time();
    /**
     * Simulation unit for time is the picosecond.
     */
    public static final Unit SIM_UNIT = Picosecond.UNIT;

    private Time() {
        super("Time", 0, 0, 1);
    }

    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.time();
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