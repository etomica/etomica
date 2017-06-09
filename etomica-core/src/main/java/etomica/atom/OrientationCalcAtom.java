/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.space3d.IOrientationFull3D;

/**
 * OrientationCalc implementation that handles a monotomic oriented molecule.
 *
 * @author Andrew Schultz
 */
public class OrientationCalcAtom implements OrientationCalc {

    public void calcOrientation(IMolecule molecule,
            IOrientationFull3D orientation) {
        orientation.E(((IAtomOriented)molecule).getOrientation());
    }

    public void setOrientation(IMolecule molecule,
            IOrientationFull3D orientation) {
        ((IAtomOriented)molecule).getOrientation().E(orientation);
    }

}
