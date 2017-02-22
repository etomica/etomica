/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.api.IAtomType;
import etomica.space.IOrientation;
import etomica.space.ISpace;
import etomica.space3d.Orientation3D;
import etomica.space3d.OrientationFull3D;

public class AtomOriented extends Atom implements
        IAtomOriented {

    public AtomOriented(ISpace space, IAtomType type) {
        this(space, type, false);
    }

    public AtomOriented(ISpace space, IAtomType type, boolean isAxisSymmetric) {
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

    private static final long serialVersionUID = 1L;
    protected final IOrientation iOrientation;
}
