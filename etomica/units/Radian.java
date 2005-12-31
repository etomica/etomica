package etomica.units;

import java.io.ObjectStreamException;

/**
 * Simulation unit for the measure of an angle.
 */
public final class Radian extends SimpleUnit {

    /**
     * Singleton instance of this unit.
     */
    public static final Radian UNIT = new Radian();

    private Radian() {
        super(Angle.DIMENSION, 1.0, "radians", "rad", Prefix.ALLOWED);
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