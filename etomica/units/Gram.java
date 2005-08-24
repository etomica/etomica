package etomica.units;
import etomica.util.Constants;

public final class Gram extends BaseUnit.Mass {

  /**
   * Convenience instance of this unit to permit unit to be assigned
   * without creating a new instance of it each time
   */
    public static final Gram UNIT = new Gram();
    
    public Gram() {
        super(
        	Constants.AVOGADRO, //conversion from grams to Daltons
        	"grams",
        	"g"
        	);   
    }
}