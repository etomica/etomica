package etomica.units;

/**
 * Simulation unit of time.
 */
public final class Picosecond extends BaseUnit.Time {

  /**
   * Singleton instance of this unit.
   */
    public static final Picosecond UNIT = new Picosecond();
  
    private Picosecond() {
        super(
        	1.0, 
        	"picoseconds",
        	"ps",
        	Prefix.NOT_ALLOWED
        	);
    }
}