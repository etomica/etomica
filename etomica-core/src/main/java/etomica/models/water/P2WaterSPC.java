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

    public P2WaterSPC(Space space) {
        super(space, 3.1670, Kelvin.UNIT.toSim(78.23), 
                Electron.UNIT.toSim(-0.82), Electron.UNIT.toSim(0.41));
    }

    private static final long serialVersionUID = 1L;
}
