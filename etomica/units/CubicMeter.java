package etomica.units;
import etomica.Constants;

public final class CubicMeter extends BaseUnit.Volume {

  /**
   * Convenience instance of this unit to permit unit to be assigned
   * without creating a new instance of it each time
   */
    public static final CubicMeter UNIT = new CubicMeter();
    
    public CubicMeter() {
        to = 1e+30; //conversion from meters^3 to Angstroms^3
        from = 1.0/to;
        name = "cubic meters";
        symbol = "m^3";   
    }
}