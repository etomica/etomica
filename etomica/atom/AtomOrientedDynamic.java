/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.api.IAtomType;
import etomica.api.IVectorMutable;
import etomica.space.IOrientation;
import etomica.space.ISpace;
import etomica.space3d.Orientation3D;
import etomica.space3d.OrientationFull3D;

public class AtomOrientedDynamic extends AtomLeafDynamic implements
        IAtomOrientedKinetic {

    private static final long serialVersionUID = 1L;
    public AtomOrientedDynamic(ISpace space, IAtomType type) {
        this(space, type, false);
    }
    
    public AtomOrientedDynamic(ISpace space, IAtomType type, boolean isAxisSymmetric) {
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
        angularVelocity = space.makeVector();  //XXX wrong! see https://rheneas.eng.buffalo.edu/bugzilla/show_bug.cgi?id=128
    }

    public IVectorMutable getAngularVelocity() {
        return angularVelocity;
    }

    public IOrientation getOrientation() {
        return iOrientation;
    }

    protected final IOrientation iOrientation;
    protected final IVectorMutable angularVelocity;
}
