
package etomica.models.water;

import etomica.space.ISpace;
import etomica.units.Calorie;
import etomica.units.Electron;
import etomica.units.Mole;

/** 
 * TIP4P potential for water.  All the real work happens in P2Water4P.
 */
public class P2WaterTIP4P extends P2Water4P {

    // U = A/r^12 + C/r^6
    private static double A = 600e3; // kcal A^12 / mol
    private static double C = 610; // kcal A^6 / mol
    private static double s6 = A/C;
    public static double s = Math.pow(s6, 1.0/6.0);
    public static double e = Mole.UNIT.fromSim(Calorie.UNIT.toSim(C/s6*1000))/4;
    
    public P2WaterTIP4P(ISpace space) {
	    super(space, s, e, Electron.UNIT.toSim(-1.04),
	            Electron.UNIT.toSim(0.52));
    }

    private static final long serialVersionUID = 1L;
}