package etomica.models.water;

import etomica.units.Calorie;
import etomica.units.Electron;
import etomica.units.Mole;

public class P2WaterTIP4P2005 {
    // U(r) = A/r^12 - C/r^6
    private static final double A = 731.3e3; // kcal Angstrom^12 / mol
    private static final double C = 736.0;   // kcal Angstrom^6 / mol
    private static final double s6 = A / C;

    // sigma_OO = 3.1589 Angstrom
    public static final double s = Math.pow(s6, 1.0 / 6.0);

    // epsilon_OO = 0.1852 kcal/mol
    public static final double e = Mole.UNIT.fromSim(Calorie.UNIT.toSim(C / s6 * 1000.0)) / 4.0;

    // Hydrogen charge = +0.5564 e
    public static final double qH = Electron.UNIT.toSim(0.5564);
}
