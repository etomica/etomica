/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.spin;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.config.Configuration;


/**
 * Configuration used to align all atoms in a spin system, so all
 * point in the same direction.  Spin direction is given by the 
 * atom's position vector.
 *
 * @author David Kofke
 *
 */
public class ConfigurationAligned implements Configuration, java.io.Serializable {

    public ConfigurationAligned() {
    }

    /**
     * Sets all spins to be aligned in the +x direction
     */
    public void initializeCoordinates(Box box) {
        IAtomList leafAtoms = box.getLeafList();
        for (int i = 0; i<leafAtoms.size(); i++) {
            Vector spin = leafAtoms.get(i).getPosition();
            spin.E(0.0);
            spin.setX(0,1.0);
        }
    }

    private static final long serialVersionUID = 2L;
}
