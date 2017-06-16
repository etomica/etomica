/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;

import etomica.space.Vector;

public interface IMoleculePositioned extends IMolecule {

    /**
     * Returns the position of the IMolecule.  Modifying the returned Vector
     * will alter the IMolecule's position, but not any positions of the
     * IMolecule's child IAtoms.
     */
    public Vector getPosition();
}
