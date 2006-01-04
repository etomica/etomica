package etomica.units;
import java.io.ObjectStreamException;

/**
 * The dyne unit of force, equal to 1 g-cm/s^2.  This
 * is the standard unit of force in the CGS unit system.
 * It is equal to 10^-5 newtons, or approximately 6.022e7 simulation force units. 
 */
public final class Dyne extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Dyne UNIT = new Dyne();
    
    private Dyne() {
        super(Force.DIMENSION,
        	Newton.UNIT.toSim(1.e-5), //6.022e7; conversion from g-cm/s^2 to Dalton-A/ps^2
        	"dynes", "dyn", Prefix.ALLOWED
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