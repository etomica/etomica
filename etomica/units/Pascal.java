package etomica.units;
import java.io.ObjectStreamException;

import etomica.util.Constants;

/**
 * The Pascal unit of pressure, equal to 1 N/m^2.
 */
public final class Pascal extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Pascal UNIT = new Pascal();
  
  /**
   * Conversion factor to/from simulation units
   */
    private Pascal() {
        super(Pressure.DIMENSION,
                Constants.AVOGADRO*1000.*1e-10*1e-24, //6.022e-8; conversion from kg/(m-s^2) to D/(A-ps^2)
                "pascals", "Pa", Prefix.ALLOWED
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