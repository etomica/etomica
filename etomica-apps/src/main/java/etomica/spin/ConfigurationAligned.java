/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.spin;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IVectorMutable;
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

    /**
     * @param space
     */
    public ConfigurationAligned() {
    }

    /**
     * Sets all spins to be aligned in the +x direction
     */
    public void initializeCoordinates(IBox box) {
        IAtomList leafAtoms = box.getLeafList();
        for (int i=0; i<leafAtoms.getAtomCount(); i++) {
            IVectorMutable spin = leafAtoms.getAtom(i).getPosition();
            spin.E(0.0);
            spin.setX(0,1.0);
        }
    }

    private static final long serialVersionUID = 2L;
}
