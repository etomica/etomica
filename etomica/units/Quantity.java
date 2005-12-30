package etomica.units;

import java.io.ObjectStreamException;

/**
 * The dimension for the quantity of discrete things (e.g. molecules). Examples
 * include Count and Mole.
 */
public final class Quantity extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Quantity();
    /**
     * Simulation unit is Count.
     */
    public static final Unit SIM_UNIT = Count.UNIT;

    private Quantity() {
        super("Quantity", 0, 0, 0, 0, 0, 1, 0);
    }

    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.quantity();
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