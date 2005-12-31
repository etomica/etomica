package etomica.units;

import java.io.ObjectStreamException;

import etomica.units.systems.UnitSystem;

/**
 * The dimension for electrical charge. Internally, charge q always represents
 * q/(4 pi epsilon0)^(1/2), so that the force via Coulomb's law is given by q1*q2/r^2. 
 * Internally, charge is a quantity having units of sqrt(force)*length, or mass^(1/2)-length^(3/2)/time, 
 * so the simulation unit for the quantity representing charge is sqrt(D-A^3)/ps.
 */
public class Charge extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Charge();
    /**
     * Simulation unit of charge.
     */
    public static final Unit SIM_UNIT = new SimpleUnit(DIMENSION, 1.0, "sim charge units", "(D-A^3/ps^2)^(1/2)/(4 pi \u03B50)^(1/2)", Prefix.NOT_ALLOWED);

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