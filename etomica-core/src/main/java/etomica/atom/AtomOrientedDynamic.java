/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Orientation3D;
import etomica.space3d.OrientationFull3D;

public class AtomOrientedDynamic extends AtomLeafDynamic implements
        IAtomOrientedKinetic {

    private static final long serialVersionUID = 1L;
    protected final IOrientation iOrientation;
    protected final Vector angularVelocity;

    public AtomOrientedDynamic(Space space, AtomType type) {
        this(space, type, false);
    }

    public AtomOrientedDynamic(Space space, AtomType type, boolean isAxisSymmetric) {
        super(space, type, space.makeVector(), space.makeVector());
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

    public AtomOrientedDynamic(Space space, AtomType type, IOrientation orientation, Vector position, Vector velocity, Vector angularVelocity) {
        super(space, type, position, velocity);
        this.iOrientation = orientation;
        this.angularVelocity = angularVelocity;
    }

    public Vector getAngularVelocity() {
        return angularVelocity;
    }

    public IOrientation getOrientation() {
        return iOrientation;
    }

    public void copyFrom(IAtom other) {
        super.copyFrom(other);
        this.iOrientation.E(((AtomOrientedDynamic) other).iOrientation);
        this.angularVelocity.E(((AtomOrientedDynamic) other).angularVelocity);
    }
}
