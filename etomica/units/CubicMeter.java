package etomica.units;

import java.io.ObjectStreamException;

public final class CubicMeter extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final CubicMeter UNIT = new CubicMeter();
    
    private CubicMeter() {
    	super(Volume.DIMENSION,
        	1e+30, //conversion from meters^3 to Angstroms^3
			"cubic meters",
        	"m^3",
        	Prefix.NOT_ALLOWED
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