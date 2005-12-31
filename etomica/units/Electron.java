package etomica.units;
import java.io.ObjectStreamException;

import etomica.util.Constants;

/**
 * Unit of charge equal to the magnitude of the charge on an electron.
 */
public final class Electron extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Electron UNIT = new Electron();
    
    private Electron() {
        super(Charge.DIMENSION, 
                1.0/Math.sqrt(4.0 * Math.PI * Constants.EPSILON_0), //372.7; conversion to (amu-A^3/ps^2)^(1/2); need only divide by sqrt(4 pi eps0) because eps0 is defined in Constants in terms of the electron charge. 
	        "elementary charges","e", Prefix.ALLOWED
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