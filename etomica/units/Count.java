package etomica.units;

/**
 * Unit class for the number, or quantity, of something.
 * This is one of the default simulation units.
 */
public final class Count extends BaseUnit.Quantity {

  /**
   * Singleton instance of this unit 
   */
    public static final Count UNIT = new Count();
    
    private Count() {
        super(
        	1.0,
        	"Count",
        	"",
        	Prefix.NOT_ALLOWED
        	);
    }
}