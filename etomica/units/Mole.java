package etomica.units;
import etomica.util.Constants;

public final class Mole extends BaseUnit.Quantity {

  /**
   * Singleton instance of this unit.
   */
    public static final Mole UNIT = new Mole();
    
    private Mole() {
        super(
        	Constants.AVOGADRO, //6.022e22; conversion from moles to count (number)
	        "mole",
        	"mol" 
        	);
    }
}