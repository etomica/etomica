/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.models.water;

import etomica.space.Space;
import etomica.units.Electron;
import etomica.units.Kelvin;

/**
 * SPC potential for water.  All the real work is done in P2Water3P.
 */
public class P2WaterSPC extends P2Water3P {

    public static final double sigmaOO = 3.1670;
    public static final double epsilonOO = Kelvin.UNIT.toSim(78.23);
    public static final double chargeO = Electron.UNIT.toSim(-0.82);
    public static final double chargeH = Electron.UNIT.toSim(0.41);

    public P2WaterSPC(Space space) {
        super(space, sigmaOO, epsilonOO, chargeO, chargeH);
    }
}
