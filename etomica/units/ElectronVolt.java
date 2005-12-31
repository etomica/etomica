package etomica.units;
import java.io.ObjectStreamException;

import etomica.util.Constants;

/**
 * The electronvolt unit of energy, equal to approximately 1.602e-19 Joules.
 */
public final class ElectronVolt extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final ElectronVolt UNIT = new ElectronVolt();

    private ElectronVolt() {
        super(Energy.DIMENSION,
         // (1.6e-19 J) (1 kg-m^2/s^2/J) (1000 g/kg) (N_avo D/g) (1e10 A/m)^2 (1e-12 s/ps)^2
        	Constants.AVOGADRO*1000.*1e20*1e-24*( 1.60217653e-19), // conversion from eV to Dalton-A^2/ps^2
        	"ElectronVolts", "eV", Prefix.ALLOWED
        	);   
    }
    
    /**
     * Required to guarantee singleton when deserializing.
     * 
     * @return the singleton UNIT
     */
    private Object readResolve() throws ObjectStreamException {
        return UNIT;
    }
    
    private static final long serialVersionUID = 1;

}