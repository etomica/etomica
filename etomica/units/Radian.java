package etomica.units;

/**
 * Simulation unit for the measure of an angle.
 */
public final class Radian extends BaseUnit.Angle {

  /**
   * Singleton instance of this unit.
   */
    public static final Radian UNIT = new Radian();
    
    private Radian() {
        super(
        	1.0, 
	        "Radians",
        	"rad"
        	);   
    }
}