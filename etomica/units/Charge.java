package etomica.units;

import java.io.ObjectStreamException;

/**
 * The dimension for electrical charge. The simulation unit of charge is the
 * charge of an electron.
 */
public class Charge extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Charge();
    public static final Unit SIM_UNIT = Electron.UNIT;

    private Charge() {
        super("Charge", 0, 0, 1, 1, 0, 0, 0);
    }

    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.charge();
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