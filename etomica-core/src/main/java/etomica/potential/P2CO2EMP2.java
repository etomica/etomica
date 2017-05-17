/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Space;
import etomica.units.Electron;
import etomica.units.Kelvin;

/**
 * EPM2 model for CO2
 */
public class P2CO2EMP2 extends P2CO2EMP {

    public P2CO2EMP2(Space space) {
        super(space, 2.757, 2.8921, 3.033, Kelvin.UNIT.toSim(28.129), Kelvin.UNIT.toSim(47.588), Kelvin.UNIT.toSim(80.507), Electron.UNIT.toSim(0.6512));
    }

    private static final long serialVersionUID = 1L;
    
}
