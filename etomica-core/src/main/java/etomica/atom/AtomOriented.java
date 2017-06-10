/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space3d.Orientation3D;
import etomica.space3d.OrientationFull3D;

public class AtomOriented extends Atom implements
        IAtomOriented {

    private static final long serialVersionUID = 1L;
    protected final IOrientation iOrientation;

    public AtomOriented(Space space, AtomType type) {
        this(space, type, false);
    }

    public AtomOriented(Space space, AtomType type, boolean isAxisSymmetric) {
        super(space, type);
        if (space.D() == 3) {
            if (isAxisSymmetric) {
                iOrientation = new Orientation3D(space);
            }
            else {
                iOrientation = new OrientationFull3D(space);
            }
        }
        else {
            iOrientation = space.makeOrientation();
        }
    }

    public IOrientation getOrientation() {
        return iOrientation;
    }
}
