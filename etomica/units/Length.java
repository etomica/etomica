package etomica.units;

import java.io.ObjectStreamException;

import etomica.units.systems.UnitSystem;

/**
 * Dimension for all units of length. Simulation unit of length is the Angstrom.
 */
public final class Length extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Length();
    public static final Unit SIM_UNIT = Angstrom.UNIT;
 
    private Length() {
        super("Length", 1, 0, 0);
    }

    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.length();
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