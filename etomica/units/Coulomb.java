package etomica.units;
import java.io.ObjectStreamException;

import etomica.util.Constants;

/**
 * The coulomb unit of electrical charge.
 */
public final class Coulomb extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Coulomb UNIT = new Coulomb();
    
    private Coulomb() {
        //note that Constants.EPSILON_0 is defined in terms of electron charge units, so must convert coulomb to e
        // then we have (1 electron-unit / 1.6e-19 C) / sqrt[4 Pi eps0]  
        super(Charge.DIMENSION,
                1.0 / (1.60217653e-19 * Math.sqrt(4 * Math.PI * Constants.EPSILON_0)), //2.326e21; conversion from Coulombs to (amu-A^3/ps^2)^(1/2)
                "coulombs", "C", Prefix.ALLOWED);   
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