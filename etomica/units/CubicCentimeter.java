package etomica.units;

import java.io.ObjectStreamException;

/**
 * The cubic centimeter unit of volume, cm^3.
 */
public final class CubicCentimeter extends SimpleUnit {

    /**
     * Singleton instance of this unit.
     */
    public static final CubicCentimeter UNIT = new CubicCentimeter();

    private CubicCentimeter() {
        super(Volume.DIMENSION, 1e+24, // conversion from cm^3 to Angstroms^3
                "cubic centimeters", "cc", Prefix.NOT_ALLOWED);
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