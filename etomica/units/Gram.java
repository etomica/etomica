package simulate.units;
import simulate.Constants;

public final class Gram extends BaseUnit.Mass {

  /**
   * Convenience instance of this unit to permit unit to be assigned
   * without creating a new instance of it each time
   */
    public static final Gram UNIT = new Gram();
    
    public Gram() {
        to = Constants.AVOGADRO; //conversion from grams to Daltons
        from = 1.0/to;
        name = "grams";
        symbol = "g";   
    }
}