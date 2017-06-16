/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.models.water;

import etomica.molecule.IMoleculePositionDefinition;
import etomica.space.Space;

/** 
 * TIP4P potential for water.  All the real work is done in P2Water4P.
 */
public class P2WaterTIP4PSoft extends P2Water4PSoft {

    public P2WaterTIP4PSoft(Space space, double rCut, IMoleculePositionDefinition positionDefinition) {
        super(space, P2WaterTIP4P.s, P2WaterTIP4P.e,
                 P2WaterTIP4P.qH,rCut,positionDefinition);
    }

    private static final long serialVersionUID = 1L;
}
