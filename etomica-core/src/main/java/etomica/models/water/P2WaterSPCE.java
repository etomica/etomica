/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.models.water;

import etomica.units.Electron;
import etomica.units.Joule;
import etomica.units.Mole;

/**
 * SPC/E potential for water.  Uses P2Water3P to do all the real work.
 *
 * @author Andrew Schultz
 */
public class P2WaterSPCE {

    private final static double B = 0.3428; //(kJ/mol)^1/12 . nm
    private final static double A = 0.37122; //(kJ/mol)^1/6 . nm
    public final static double EPSILON = Mole.UNIT.fromSim(Joule.UNIT.toSim(Math.pow((A/B),12)/4))*1000;
    public final static double SIGMA = 10*B*B/A;
    public final static double QH = Electron.UNIT.toSim(0.4238);
    public final static double QO = -2 * QH;
}
