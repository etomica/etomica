package etomica.units;
import etomica.Constants;

public final class Kelvin extends BaseUnit.Temperature {

  /**
   * Convenience instance of this unit to permit unit to be assigned
   * without creating a new instance of it each time
   */
    public static final Kelvin UNIT = new Kelvin();
  
    public Kelvin() {
        to = Constants.BOLTZMANN_K; 
        from = 1.0/to;
        name = "Kelvins";
        symbol = "K";   
    }
}