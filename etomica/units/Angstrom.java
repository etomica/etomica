package etomica.units;

public final class Angstrom extends BaseUnit.Length {

  /**
   * Singleton instance of this unit
   */
    public static final Angstrom UNIT = new Angstrom();
  
    private Angstrom() {
    	super(
    		1.0,//conversion to simulation units
	        "Angstroms",
	        "\u00c5",   //unicode for the Angstrom symbol
	        Prefix.NOT_ALLOWED
	        );
    }
}