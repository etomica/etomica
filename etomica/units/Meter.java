package etomica.units;

public final class Meter extends BaseUnit.Length {

  /**
   * Singleton instance of this unit.
   */
    public static final Meter UNIT = new Meter();
    
    private Meter() {
    	super(
        	1e+10, //conversion from meters to Angstroms
	        "meters",
	        "m"
    	);
    }
}