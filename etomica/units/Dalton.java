package simulate.units;

public final class Dalton extends BaseUnit.Mass {

  /**
   * Convenience instance of this unit to permit unit to be assigned
   * without creating a new instance of it each time
   */
    public static final Dalton UNIT = new Dalton();
  
    public Dalton() {
        to = 1.0; 
        from = 1.0/to;
        name = "Daltons";
        symbol = "D";
    }
}