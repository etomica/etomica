/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.space.Space;
import etomica.space.Vector;
import etomica.spaceNd.VectorND;

public class AtomOrientedQuaternion extends Atom {

    protected final VectorND quaternion;

    public AtomOrientedQuaternion(Space space, AtomType type) {
        super(space, type);
        quaternion = new VectorND(4);
        quaternion.setX(0, 1);
    }

    public Vector getQuaternion() {
        return quaternion;
    }
}
