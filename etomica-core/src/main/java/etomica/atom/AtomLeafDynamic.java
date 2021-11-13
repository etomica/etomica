/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.space.Space;
import etomica.space.Vector;

public class AtomLeafDynamic extends Atom implements IAtomKinetic {

    protected final Vector velocity;

    public AtomLeafDynamic(Space space, AtomType type) {
        super(space, type);
        velocity = space.makeVector();
    }

    public Vector getVelocity() {
        return velocity;
    }

    public void copyCoordinatesFrom(IAtom atom) {
        super.copyCoordinatesFrom(atom);
        velocity.E(((IAtomKinetic) atom).getVelocity());
    }
}
