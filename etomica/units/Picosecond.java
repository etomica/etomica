package etomica.units;

public final class Picosecond extends BaseUnit.Time {

  /**
   * Convenience instance of this unit to permit unit to be assigned
   * without creating a new instance of it each time
   */
    public static final Picosecond UNIT = new Picosecond();
  
    public Picosecond() {
        prefixAllowed = false;
        to = 1.0; 
        from = 1.0/to;
        name = "picoseconds";
        symbol = "ps";
    }
}