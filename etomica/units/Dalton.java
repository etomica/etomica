package etomica.units;

import java.io.ObjectStreamException;

public final class Dalton extends SimpleUnit {

    /**
     * Singleton instance of this unit.
     */
    public static final Dalton UNIT = new Dalton();

    private Dalton() {
        super(Mass.DIMENSION, 1.0, "daltons", "Da", Prefix.ALLOWED);
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