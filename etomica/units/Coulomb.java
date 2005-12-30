package etomica.units;
import java.io.ObjectStreamException;

import etomica.util.Constants;

public final class Coulomb extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Coulomb UNIT = new Coulomb();
    
    private Coulomb() {
        //should check this more carefully
        // sqrt[ 1/(4 Pi epsilon0) * (A/m)^2 * (g/kg) * (D/g) * (A/m) * (s/ps)^2 ]
        super(Charge.DIMENSION,
                Math.sqrt(1/4/Math.PI/8.854e-12*1e20*1000*Constants.AVOGADRO*1e10*1e-24), //2.326e23; conversion from Coulombs to (amu-A^3/ps^2)^(1/2)
                "Coulombs", "C", Prefix.ALLOWED);   
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