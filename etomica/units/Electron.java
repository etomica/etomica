package etomica.units;
import etomica.util.Constants;

/**
 * Unit of charge equal to the magnitude of the charge on an electron
 */
public final class Electron extends BaseUnit.Charge {

  /**
   * Singleton instance of this unit.
   */
    public static final Electron UNIT = new Electron();
    
    private Electron() {
        super(
        	4.803e-10*Math.sqrt(Constants.AVOGADRO*1e24*1e-24), //372.7; conversion from (electron/esu)*(g-cm^3/s^2)^(1/2) to (amu-A^3/ps^2)^(1/2)
	        "electron-charge units",
        	"e"
        	);   
    }
}