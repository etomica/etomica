package etomica.units;

import java.io.ObjectStreamException;

import etomica.util.Constants;

/**
 * The Joule unit of energy, equal to 1 N-m or 1 kg-m^2/s^2.
 */
public final class Calorie extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Calorie UNIT = new Calorie();
    
    private Calorie() {
        super(Energy.DIMENSION,
        	Constants.AVOGADRO*1000.*1e20*1e-24/0.238845896628, //6.022e22; conversion from calories to Joules kg-m^2/s^2 to Dalton-A^2/ps^2
        	"calories", "cal", Prefix.ALLOWED
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
