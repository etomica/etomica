package etomica.units;

public final class CubicMeter extends BaseUnit.Volume {

  /**
   * Singleton instance of this unit.
   */
    public static final CubicMeter UNIT = new CubicMeter();
    
    private CubicMeter() {
    	super(
        	1e+30, //conversion from meters^3 to Angstroms^3
			"cubic meters",
        	"m^3",
        	Prefix.NOT_ALLOWED
        	);   
    }
}