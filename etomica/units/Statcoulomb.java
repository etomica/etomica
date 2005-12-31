package etomica.units;
import java.io.ObjectStreamException;

import etomica.util.Constants;

/**
 * The statcoulomb unit of electrical charge, which is the 
 * standard unit of charge in the CGS unit system.  Also known as
 * the electrostatic unit (esu), or the franklin.  Equal to
 * approximately 3.3356e-10 coulombs, or sqrt(N_avogadro) in
 * simulation units.
 *
 */
public final class Statcoulomb extends SimpleUnit {

  /**
   * Convenience instance of this unit to permit unit to be assigned
   * without creating a new instance of it each time
   */
    public static final Statcoulomb UNIT = new Statcoulomb();
    
    private Statcoulomb() {
        super(Charge.DIMENSION,
                Math.sqrt(Constants.AVOGADRO*1e24*1e-24), //7.76e11; conversion from (g-cm^3/s^2)^(1/2) to (D-A^3/ps^2)^(1/2)
                "electrostatic units", "esu", Prefix.ALLOWED
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