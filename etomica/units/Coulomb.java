package simulate.units;
import simulate.Constants;

public final class Coulomb extends BaseUnit.Charge {

  /**
   * Convenience instance of this unit to permit unit to be assigned
   * without creating a new instance of it each time
   */
    public static final Coulomb UNIT = new Coulomb();
    
    public Coulomb() {
        //should check this more carefully
        // sqrt[ 1/(4 Pi epsilon0) * (A/m)^2 * (g/kg) * (D/g) * (A/m) * (s/ps)^2 ]
        to = Math.sqrt(1/4/Math.PI/8.854e-12*1e20*1000*Constants.AVOGADRO*1e10*1e-24); //2.326e23; conversion from Coulombs to (amu-A^3/ps^2)^(1/2)
        from = 1.0/to;
        name = "Coulombs";
        symbol = "C";   
    }
}