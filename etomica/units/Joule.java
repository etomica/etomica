package etomica.units;
import java.io.ObjectStreamException;

import etomica.util.Constants;

public final class Joule extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Joule UNIT = new Joule();
    
    private Joule() {
        super(Energy.DIMENSION,
        	Constants.AVOGADRO*1000.*1e20*1e-24, //6.022e22; conversion from kg-m^2/s^2 to Dalton-A^2/ps^2
        	"Joules", "J", Prefix.ALLOWED
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