package etomica.units;
import etomica.Constants;

public final class Radian extends BaseUnit.Angle {

  /**
   * Convenience instance of this unit to permit unit to be assigned
   * without creating a new instance of it each time
   */
    public static final Radian UNIT = new Radian();
    
    public Radian() {
        to = 1.0; 
        from = 1.0/to;
        name = "Radians";
        symbol = "rad";   
    }
}