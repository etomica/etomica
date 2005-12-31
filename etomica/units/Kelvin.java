package etomica.units;
import java.io.ObjectStreamException;

import etomica.util.Constants;

public final class Kelvin extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Kelvin UNIT = new Kelvin();
  
    private Kelvin() {
        super(Temperature.DIMENSION, 
        	Constants.BOLTZMANN_K,//convert to simulation energy units by multiplying by Boltzmann's constant
        	"kelvins", "K", Prefix.ALLOWED
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