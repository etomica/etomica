package etomica.units;
import java.io.ObjectStreamException;

import etomica.util.Constants;

/**
 * Unit of charge equal to the magnitude of the charge on an electron
 */
public final class Electron extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Electron UNIT = new Electron();
    
    private Electron() {
        super(Charge.DIMENSION, 
        	4.803e-10*Math.sqrt(Constants.AVOGADRO*1e24*1e-24), //372.7; conversion from (electron/esu)*(g-cm^3/s^2)^(1/2) to (amu-A^3/ps^2)^(1/2)
	        "electron-charge units","e", Prefix.ALLOWED
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