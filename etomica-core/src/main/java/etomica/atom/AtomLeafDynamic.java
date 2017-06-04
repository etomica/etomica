/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.space.Vector;
import etomica.space.Space;

public class AtomLeafDynamic extends Atom implements IAtomKinetic {

    public AtomLeafDynamic(Space space, IAtomType type) {
        super(space, type);
        velocity = space.makeVector();
    }
    
    public Vector getVelocity() {
        return velocity;
    }
    
    private static final long serialVersionUID = 1L;
    protected final Vector velocity;
}
