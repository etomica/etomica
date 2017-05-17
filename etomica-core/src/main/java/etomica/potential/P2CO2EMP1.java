/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Space;
import etomica.units.Electron;
import etomica.units.Kelvin;

/**
 * EPM model for CO2.  we'll call it EPM1 since we used EPM for the superclass.
 */
public class P2CO2EMP1 extends P2CO2EMP {

    public P2CO2EMP1(Space space) {
        super(space, 2.785, 2.921, 3.064, Kelvin.UNIT.toSim(28.999), Kelvin.UNIT.toSim(49.060), Kelvin.UNIT.toSim(82.9997), Electron.UNIT.toSim(0.6645));
    }

    private static final long serialVersionUID = 1L;
    
}
