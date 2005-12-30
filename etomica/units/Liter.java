package etomica.units;

import java.io.ObjectStreamException;

public final class Liter extends SimpleUnit {

  /**
   * Convenience instance of this unit to permit unit to be assigned
   * without creating a new instance of it each time
   */
    public static final Liter UNIT = new Liter();
    
    public Liter() {
        super(Volume.DIMENSION,
        	1e+27, //conversion from liters to Angstroms^3 (10^-3 m^3/liter * (10^10 Angstroms/meter)^3)
        	"liters", "l", Prefix.ALLOWED
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