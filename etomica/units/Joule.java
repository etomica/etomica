package simulate.units;
import simulate.Constants;

public final class Joule extends BaseUnit.Energy {

  /**
   * Convenience instance of this unit to permit unit to be assigned
   * without creating a new instance of it each time
   */
    public static final Joule UNIT = new Joule();
    
    public Joule() {
        to = Constants.AVOGADRO*1000.*1e20*1e-24; //6.022e22; conversion from kg-m^2/s^2 to Dalton-A^2/ps^2
        from = 1.0/to;
        name = "Joules";
        symbol = "J";   
    }
}