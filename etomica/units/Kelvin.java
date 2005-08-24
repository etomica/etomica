package etomica.units;
import etomica.util.Constants;

public final class Kelvin extends BaseUnit.Temperature {

  /**
   * Singleton instance of this unit.
   */
    public static final Kelvin UNIT = new Kelvin();
  
    private Kelvin() {
        super(
        	Constants.BOLTZMANN_K,//convert to simulation energy units by multiplying by Boltzmann's constant
        	"Kelvins",
        	"K"
        	);  
    }
}