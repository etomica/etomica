package etomica.units;

import java.io.ObjectStreamException;

/**
 * Base for all units of time.  
 * Simulation unit for time is the picosecond.
 */
public final class Time extends Dimension {

    public static final Dimension DIMENSION = new Time();
    public static final Unit SIM_UNIT = Picosecond.UNIT;
     
    private Time() {
        super("Time", 0,0,1);
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