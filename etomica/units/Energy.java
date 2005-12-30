package etomica.units;

import java.io.ObjectStreamException;

import etomica.units.systems.UnitSystem;

/**
 * Base for all energy units. Simulation unit of energy is D-A^2/ps^2
 */
public final class Energy extends Dimension {

    public static final Dimension DIMENSION = new Energy();
    public static final Unit SIM_UNIT = new SimpleUnit(DIMENSION, 1.0, "sim energy units", "D-A^2/ps^2", Prefix.NOT_ALLOWED);

    private Energy() {
        super("Energy", 2, 1, -2);
    }
    
    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.energy();
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