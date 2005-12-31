package etomica.units;
import java.io.ObjectStreamException;

import etomica.util.Constants;

/**
 * The Newton unit of force, equal to 1 kg-m/s^2.
 */
public final class Newton extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Newton UNIT = new Newton();
    
    private Newton() {
        super(Force.DIMENSION,
        	Constants.AVOGADRO*1000.*1e10*1e-24, //6.022e12; conversion from kg-m/s^2 to Dalton-A/ps^2
        	"newtons", "N", Prefix.ALLOWED
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