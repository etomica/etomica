package etomica.units;

import java.io.ObjectStreamException;

/**
 * Base unit for electrical dipole moment. Simulation unit is electron-Angstrom.
 */
public class Dipole extends Dimension {

    public static final Dimension DIMENSION = new Dipole();
    public static final Unit SIM_UNIT = new SimpleUnit(DIMENSION, 1, "sim dipole units", "e-\u00c5", Prefix.NOT_ALLOWED);

    private Dipole() {
        super("Dipole", 1, 0, 1, 1, 0, 0, 0);
    }
    
    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.dipole();
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