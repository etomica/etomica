package etomica.units;
import java.io.ObjectStreamException;

import etomica.util.Constants;

/**
 * The bar unit of pressure, equal to 10^5 N/m^2.
 */
public final class Bar extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Bar UNIT = new Bar();
  
  /**
   * Conversion factor to/from simulation units
   */
    private Bar() {
        super(Pressure.DIMENSION,
                1e5*Constants.AVOGADRO*1000.*1e-10*1e-24, //6.022e-3; conversion from 10^5 kg/(m-s^2) to amu/(A-ps^2)
                "bars", "bar", Prefix.ALLOWED
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