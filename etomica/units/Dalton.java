package etomica.units;

public final class Dalton extends BaseUnit.Mass {

  /**
   * Singleton instance of this unit.
   */
    public static final Dalton UNIT = new Dalton();
  
    private Dalton() {
        super(
			1.0,
			"Daltons",
        	"D"
        	);
    }
}