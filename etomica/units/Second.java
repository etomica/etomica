package etomica.units;

public final class Second extends BaseUnit.Time {

  /**
   * Convenience instance of this unit to permit unit to be assigned
   * without creating a new instance of it each time
   */
    public static final Second UNIT = new Second();
    
    public Second() {
        to = 1e+12; //conversion from seconds to picoseconds
        from = 1.0/to;
        name = "seconds";
        symbol = "s";   
    }
}