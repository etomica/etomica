/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.space.Space;

/**
 * Superclass for all Potential classes, which define how the atoms in the
 * system interact with each other.
 *
 * @author David Kofke
 */
 
public abstract class Potential implements IPotentialAtomic {
    
	protected final Space space;

    /**
     * General constructor for a potential instance
     */
    public Potential(Space space) {
        this.space = space;
    }

    public abstract double getRange();
    
    /**
     * Returns the interaction energy between the given atoms.  There might be
     * 0, 1, 2 or more atoms in the AtomSet.
     */
    public abstract double energy(IAtomList atoms);

}
