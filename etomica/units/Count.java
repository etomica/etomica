package simulate.units;

/**
 * Unit class for the number, or quantity, of something.
 * This is one of the default simulation units.
 */
public final class Count extends BaseUnit.Quantity {

  /**
   * Convenience instance of this unit to permit unit to be assigned
   * without creating a new instance of it each time
   */
    public static final Count UNIT = new Count();
    
    public Count() {
        to = 1.0;
        from = 1.0/to;
        name = "Count";
        symbol = ""; 
        prefixAllowed = false;
    }
}