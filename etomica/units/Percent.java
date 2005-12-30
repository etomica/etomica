package etomica.units;

import java.io.ObjectStreamException;

/**
 * Decimal representation of something that represents the fractional 
 * amount of a whole (e.g., mole fraction) as a percentage value typically
 * between 0 and 100.
 */
public final class Percent extends SimpleUnit {

  /**
   * Singleton instance of this unit 
   */
	public static final Percent UNIT = new Percent();

	private Percent() {
       super(Fraction.DIMENSION,
        	0.01,
        	"Percent", "%", Prefix.NOT_ALLOWED
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
