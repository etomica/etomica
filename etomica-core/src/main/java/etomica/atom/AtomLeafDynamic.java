/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.space.Space;
import etomica.space.Vector;

public class AtomLeafDynamic extends Atom implements IAtomKinetic {

    private static final long serialVersionUID = 1L;
    protected final Vector velocity;

    public AtomLeafDynamic(Space space, AtomType type, Vector position, Vector velocity) {
        super(space, type, position);
        this.velocity = velocity;
    }

    public Vector getVelocity() {
        return velocity;
    }

    @Override
    public void copyFrom(IAtom atom) {
        super.copyFrom(atom);
        this.velocity.E(((AtomLeafDynamic) atom).getVelocity());
    }
}
