package etomica.units;

public final class Degree extends BaseUnit.Angle {

  /**
   * Single instance of this unit.
   */
    public static final Degree UNIT = new Degree();
    
    private Degree() {
        super(
        	Math.PI/180., //conversion from degrees to radians
        	"degrees",
        	"\u00B0" //unicode for the degree symbol
        	);
    }
}