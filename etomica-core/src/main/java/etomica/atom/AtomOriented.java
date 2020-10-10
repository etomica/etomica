/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Orientation3D;
import etomica.space3d.OrientationFull3D;

public class AtomOriented extends Atom implements
        IAtomOriented {

    private static final long serialVersionUID = 1L;
    protected final IOrientation iOrientation;

    public AtomOriented(Space space, AtomType type) {
        this(space, type, space.makeVector(), space.makeOrientation());
    }

    public AtomOriented(Space space, AtomType type, Vector position, IOrientation orientation) {
        super(space, type, position);
        iOrientation = orientation;
    }

    public IOrientation getOrientation() {
        return iOrientation;
    }

    @Override
    public void copyFrom(IAtom atom) {
        super.copyFrom(atom);
        this.iOrientation.E(((AtomOriented) atom).iOrientation);
    }
}
