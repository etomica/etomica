package simulate.units;
import simulate.Constants;

public final class Degree extends BaseUnit.Angle {

  /**
   * Convenience instance of this unit to permit unit to be assigned
   * without creating a new instance of it each time
   */
    public static final Degree UNIT = new Degree();
    
    public Degree() {
        to = Math.PI/180.; //conversion from degrees to radians
        from = 1.0/to;
        name = "degrees";
        symbol = "\u00B0"; //unicode for the degree symbol
    }
}