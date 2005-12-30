package etomica.units;
import java.io.ObjectStreamException;

import etomica.util.Constants;

public final class ElectronVolt extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final ElectronVolt UNIT = new ElectronVolt();
    
    private ElectronVolt() {
        super(Energy.DIMENSION,
        	Constants.AVOGADRO*1000.*1e20*1e-24*(1.6e-19), // conversion from eV to Dalton-A^2/ps^2
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