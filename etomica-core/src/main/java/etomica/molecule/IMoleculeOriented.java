/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;

import etomica.space.IOrientation;

/**
 * Interface for a IAtom that includes an IVector that defines the atom's
 * position.
 */
public interface IMoleculeOriented extends IMoleculePositioned {
    
    /**
     * Returns the orientation of the IAtom.  Modifying the returned IVector will
     * alter the IAtom's position.
     */
    public IOrientation getOrientation();

}
