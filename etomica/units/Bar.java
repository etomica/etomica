package etomica.units;
import etomica.Constants;

public final class Bar extends BaseUnit.Pressure implements BaseUnit.D3 {

  /**
   * Convenience instance of this unit to permit unit to be assigned
   * without creating a new instance of it each time
   */
    public static final Bar UNIT = new Bar();
  
  /**
   * Conversion factor to/from simulation units
   */
    public Bar() {
        to = 1e5*Constants.AVOGADRO*1000.*1e-10*1e-24; //6.022e-3; conversion from 10^5 kg/(m-s^2) to amu/(A-ps^2)
        from = 1.0/to;
        name = "bars";
        symbol = "bar";   
    }
}