/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomArrayList;
import etomica.molecule.IMoleculeList;

/**
 * Intermolecular potential that wraps an atomic potential.  The potential is
 * assumed to be monatomic.
 * 
 * @author Andrew Schultz
 */
public class PotentialMolecularMonatomic implements IPotentialMolecular {
    
    protected final IPotentialAtomic potentialAtomic;
    protected final AtomArrayList atoms;

    public PotentialMolecularMonatomic(IPotentialAtomic potentialAtomic, int nBody) {
        super();
        this.potentialAtomic = potentialAtomic;
        atoms = new AtomArrayList(nBody);
    }

    public double energy(IMoleculeList molecules) {
        atoms.clear();
        for (int i=0; i<atoms.size(); i++) {
            atoms.add(molecules.get(i).getChildList().get(0));
        }
        return potentialAtomic.energy(atoms);
    }
}
