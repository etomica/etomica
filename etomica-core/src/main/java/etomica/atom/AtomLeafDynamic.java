/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.api.IAtomKinetic;
import etomica.api.IAtomType;
import etomica.api.IVectorMutable;
import etomica.space.ISpace;

public class AtomLeafDynamic extends Atom implements IAtomKinetic {

    public AtomLeafDynamic(ISpace space, IAtomType type) {
        super(space, type);
        velocity = space.makeVector();
    }
    
    public IVectorMutable getVelocity() {
        return velocity;
    }
    
    private static final long serialVersionUID = 1L;
    protected final IVectorMutable velocity;
}
