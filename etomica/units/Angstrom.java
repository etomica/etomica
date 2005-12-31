package etomica.units;

import java.io.ObjectStreamException;

public final class Angstrom extends SimpleUnit {

    /**
     * Singleton instance of this unit
     */
    public static final Angstrom UNIT = new Angstrom();

    private Angstrom() {
        super(Length.DIMENSION, 1.0,// conversion to simulation units
                "angstroms", "\u00c5", // unicode for the Angstrom symbol
                Prefix.NOT_ALLOWED);
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