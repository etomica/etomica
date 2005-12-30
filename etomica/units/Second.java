package etomica.units;

import java.io.ObjectStreamException;

public final class Second extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Second UNIT = new Second();
    
    private Second() {
        super(Time.DIMENSION,
                1e+12, //conversion from seconds to picoseconds
                "seconds", "s", Prefix.ALLOWED
        	);   
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