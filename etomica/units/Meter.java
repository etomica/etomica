package etomica.units;

import java.io.ObjectStreamException;

/**
 * The meter unit of length, equal to 10^10 angstroms.
 */
public final class Meter extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Meter UNIT = new Meter();
    
    private Meter() {
    	    super(Length.DIMENSION,
    	            1e+10, //conversion from meters to Angstroms
    	            "meters", "m", Prefix.ALLOWED
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