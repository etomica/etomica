package etomica.units;

public final class Liter extends BaseUnit.Volume {

  /**
   * Convenience instance of this unit to permit unit to be assigned
   * without creating a new instance of it each time
   */
    public static final Liter UNIT = new Liter();
    
    public Liter() {
        super(
        	1e+27, //conversion from liters to Angstroms^3
        	"liters",
        	"l"
        	);   
    }
}