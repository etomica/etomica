package etomica.units;

public final class Angstrom extends BaseUnit.Length {

  /**
   * Convenience instance of this unit to permit unit to be assigned
   * without creating a new instance of it each time
   */
    public static final Angstrom UNIT = new Angstrom();
  
    public Angstrom() {
        prefixAllowed = false;
        to = 1.0; 
        from = 1.0/to;
        name = "Angstroms";
        symbol = "\u00c5";   //unicode for the Angstrom symbol
    }
}