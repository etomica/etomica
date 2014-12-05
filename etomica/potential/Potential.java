/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IPotentialAtomic;
import etomica.space.ISpace;

/**
 * Superclass for all Potential classes, which define how the atoms in the
 * system interact with each other.
 *
 * @author David Kofke
 */
 
public abstract class Potential implements IPotentialAtomic {
    
	private final int nBody;
	protected final ISpace space;

    /**
     * General constructor for a potential instance
     * @param nBody number of atoms to which this potential applies at a time;
     * for example with a pair potential nBody = 2; for a single-body potential,
     * nBody = 1.
     */
    public Potential(int nBody, ISpace space) {
        this.nBody = nBody;
        this.space = space;
    }

    public abstract double getRange();
    
    /**
     * Returns the interaction energy between the given atoms.  There might be
     * 0, 1, 2 or more atoms in the AtomSet.
     */
    public abstract double energy(IAtomList atoms);
    
    /**
     * Informs the potential of the box on which it acts. Typically this
     * requires at least that it update the nearestImageTransformer of its
     * coordinatePair (if it uses one), e.g.:
     * cPair.setNearestImageTransformer(box.boundary());
     */
    public abstract void setBox(IBox box);
    
    /**
     * The number of atoms on which the potential depends.
     */
    public final int nBody() {return nBody;}
    

}