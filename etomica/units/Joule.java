package etomica.units;
import etomica.Constants;

public final class Joule extends BaseUnit.Energy {

  /**
   * Singleton instance of this unit.
   */
    public static final Joule UNIT = new Joule();
    
    private Joule() {
        super(
        	Constants.AVOGADRO*1000.*1e20*1e-24, //6.022e22; conversion from kg-m^2/s^2 to Dalton-A^2/ps^2
        	"Joules",
        	"J"
        	);   
    }
}