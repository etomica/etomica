package etomica.units;
import etomica.Constants;

public final class Liter extends BaseUnit.Volume {

  /**
   * Convenience instance of this unit to permit unit to be assigned
   * without creating a new instance of it each time
   */
    public static final Liter UNIT = new Liter();
    
    public Liter() {
        to = 1e+27; //conversion from liters to Angstroms^3
        from = 1.0/to;
        name = "liters";
        symbol = "l";   
    }
}