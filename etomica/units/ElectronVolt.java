package etomica.units;
import etomica.util.Constants;

public final class ElectronVolt extends BaseUnit.Energy {

  /**
   * Singleton instance of this unit.
   */
    public static final ElectronVolt UNIT = new ElectronVolt();
    
    private ElectronVolt() {
        super(
        	Constants.AVOGADRO*1000.*1e20*1e-24*(1.6e-19), // conversion from eV to Dalton-A^2/ps^2
        	"ElectronVolts",
        	"eV"
        	);   
    }
}