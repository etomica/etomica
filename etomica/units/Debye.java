package etomica.units;
import java.io.ObjectStreamException;

import etomica.util.Constants;

public final class Debye extends SimpleUnit {

  /**
   * Single instance of this unit.
   */
    public static final Debye UNIT = new Debye();
    
    private Debye() {
        super(Dipole.DIMENSION,
                Math.sqrt(Constants.AVOGADRO*1e40*1e-24)*1e-18, //77.6; conversion from (g-cm^5/s^2)^(1/2) to 10^-18 * (Daltons-A^5/ps^2)^(1/2) (Debye)
        	        "Debyes",
        	        "D",
                 Prefix.ALLOWED);   
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