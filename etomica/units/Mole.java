package etomica.units;
import etomica.Constants;

public final class Mole extends BaseUnit.Quantity {

  /**
   * Convenience instance of this unit to permit unit to be assigned
   * without creating a new instance of it each time
   */
    public static final Mole UNIT = new Mole();
    
    public Mole() {
        to = Constants.AVOGADRO; //6.022e22; conversion from moles to count (number)
        from = 1.0/to;
        name = "mole";
        symbol = "mol"; 
        prefixAllowed = true;
    }
}