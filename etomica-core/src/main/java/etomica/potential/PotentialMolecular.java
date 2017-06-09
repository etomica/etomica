/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.box.Box;
import etomica.atom.IMoleculeList;
import etomica.api.IPotentialMolecular;
import etomica.space.Space;

/**
 * Superclass for all Potential classes, which define how the atoms in the
 * system interact with each other.
 *
 * @author David Kofke
 */
 
public abstract class PotentialMolecular implements IPotentialMolecular {
    
	protected final int nBody;
	protected final Space space;

    /**
     * General constructor for a potential instance
     * @param nBody number of atoms to which this potential applies at a time;
     * for example with a pair potential nBody = 2; for a single-body potential,
     * nBody = 1.
     */
    public PotentialMolecular(int nBody, Space space) {
        this.nBody = nBody;
        this.space = space;
    }

    public abstract double getRange();
    
    /**
     * Returns the interaction energy between the given molecules.  There might
     * be 0, 1, 2 or more atoms in the IMoleculeList.
     */
    public abstract double energy(IMoleculeList molecules);
    
    /**
     * Informs the potential of the box on which it acts. Typically this
     * requires at least that it update the nearestImageTransformer of its
     * coordinatePair (if it uses one), e.g.:
     * cPair.setNearestImageTransformer(box.boundary());
     */
    public abstract void setBox(Box box);
    
    /**
     * The number of atoms on which the potential depends.
     */
    public final int nBody() {return nBody;}
    

}
