package etomica.units;
import etomica.Constants;

public final class Debye extends BaseUnit.Dipole {

  /**
   * Convenience instance of this unit to permit unit to be assigned
   * without creating a new instance of it each time
   */
    public static final Debye UNIT = new Debye();
    
    public Debye() {
        to = Math.sqrt(Constants.AVOGADRO*1e40*1e-24)*1e-18; //77.6; conversion from (g-cm^5/s^2)^(1/2) to 10^-18 * (Daltons-A^5/ps^2)^(1/2) (Debye)
        from = 1.0/to;
        name = "Debyes";
        symbol = "D";   
    }
}