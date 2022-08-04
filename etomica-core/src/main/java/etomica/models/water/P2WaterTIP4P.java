/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.models.water;

import etomica.units.Calorie;
import etomica.units.Electron;
import etomica.units.Mole;

/** 
 * TIP4P potential for water.  All the real work happens in P2Water4P.
 */
public class P2WaterTIP4P {

    // U = A/r^12 + C/r^6
    private static final double A = 600e3; // kcal A^12 / mol
    private static final double C = 610; // kcal A^6 / mol
    private static final double s6 = A/C;
    public static final double s = Math.pow(s6, 1.0/6.0);
    public static final double e = Mole.UNIT.fromSim(Calorie.UNIT.toSim(C/s6*1000))/4;
    public static final double qH = Electron.UNIT.toSim(0.52);
}
