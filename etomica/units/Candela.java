package etomica.units;

import java.io.ObjectStreamException;

/**
 * The candela is the luminous intensity, in a given direction, of a source that
 * emits monochromatic radiation of frequency 540 x 1012 hertz and that has a
 * radiant intensity in that direction of 1/683 watt per steradian.
 */
public final class Candela extends SimpleUnit {

    /**
     * Singleton instance of this unit.
     */
    public static final Candela UNIT = new Candela();

    private Candela() {
        super(LuminousIntensity.DIMENSION, 1.0, "Candela", "cd", Prefix.ALLOWED);
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