package etomica.units;

import java.io.ObjectStreamException;

import etomica.units.systems.UnitSystem;

/**
 * Simulation unit of (2D) pressure is (D-A/ps^2)/A = D/ps^2
 */
public final class Pressure2D extends Dimension {

    public static final Dimension DIMENSION = new Pressure2D();
    public static final Unit SIM_UNIT = new SimpleUnit(DIMENSION, 1, "sim 2-D pressure units", "D/ps^2", Prefix.NOT_ALLOWED);

    private Pressure2D() {
        super("2-D Pressure", 0, 1, -2, 0, 0, 0, 0);
    }
    
    public Unit getUnit(UnitSystem unitSystem) {
        Unit[] units = new Unit[] {Pressure.DIMENSION.getUnit(unitSystem), Length.DIMENSION.getUnit(unitSystem)};
        double[] exponents = new double[] {1.0, 1.0};
        return new CompoundUnit(units, exponents);
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