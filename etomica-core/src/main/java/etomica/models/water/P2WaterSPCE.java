/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.models.water;

import etomica.space.Space;
import etomica.units.Electron;
import etomica.units.Kelvin;

/** 
 * SPC/E potential for water.  Uses P2Water3P to do all the real work.
 * 
 * @author Andrew Schultz
 */
public class P2WaterSPCE extends P2Water3P {

    public P2WaterSPCE(Space space) {
	    super(space, 3.1670, Kelvin.UNIT.toSim(78.21), 
	            Electron.UNIT.toSim(-0.8476), Electron.UNIT.toSim(0.4238));
    }
    
    private static final long serialVersionUID = 1L;
}
