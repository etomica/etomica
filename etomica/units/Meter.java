package etomica.units;
import etomica.Constants;

public final class Meter extends BaseUnit.Length {

  /**
   * Convenience instance of this unit to permit unit to be assigned
   * without creating a new instance of it each time
   */
    public static final Meter UNIT = new Meter();
    
    public Meter() {
        to = 1e+10; //conversion from meters to Angstroms
        from = 1.0/to;
        name = "meters";
        symbol = "m";   
    }
}