/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.models.water;

import etomica.atom.IAtomPositionDefinition;
import etomica.space.ISpace;
import etomica.units.Electron;
import etomica.units.Kelvin;

/** 
 * TIP4P potential for water.  All the real work is done in P2Water4P.
 */
public class P2WaterTIP4PSoft extends P2Water4PSoft {

    public P2WaterTIP4PSoft(ISpace space,double rCut,IAtomPositionDefinition positionDefinition) {
        super(space, 3.1540, Kelvin.UNIT.toSim(78.02), 
                Electron.UNIT.toSim(-1.04), Electron.UNIT.toSim(0.52),rCut,positionDefinition);
    }

    private static final long serialVersionUID = 1L;
}
