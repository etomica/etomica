package etomica.units;
import etomica.Constants;

public final class Bar extends BaseUnit.Pressure implements BaseUnit.D3 {

  /**
   * Singleton instance of this unit.
   */
    public static final Bar UNIT = new Bar();
  
  /**
   * Conversion factor to/from simulation units
   */
    private Bar() {
        super(
        	1e5*Constants.AVOGADRO*1000.*1e-10*1e-24, //6.022e-3; conversion from 10^5 kg/(m-s^2) to amu/(A-ps^2)
        	"bars",
        	"bar"
        	);   
    }
}