package etomica.units;
import java.io.ObjectStreamException;

import etomica.util.Constants;

/**
 * The mole unit of quantity, approximately equal to 6.022e23 simulation units.
 */
public final class Mole extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Mole UNIT = new Mole();
    
    private Mole() {
        super(Quantity.DIMENSION,
               Constants.AVOGADRO, //6.022e23; conversion from moles to count (number)
               "moles",
               "mol",
               Prefix.ALLOWED
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