/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.api.IPotentialAtomic;
import etomica.atom.AtomArrayList;
import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;

/**
 * Intermolecular potential that wraps an atomic potential.  The potential is
 * assumed to be monatomic.
 * 
 * @author Andrew Schultz
 */
public class PotentialMolecularMonatomic extends PotentialMolecular {
    
    protected final IPotentialAtomic potentialAtomic;
    protected final AtomArrayList atoms;

    public PotentialMolecularMonatomic(Space space, IPotentialAtomic potentialAtomic) {
        super(potentialAtomic.nBody(), space);
        this.potentialAtomic = potentialAtomic;
        atoms = new AtomArrayList(nBody);
    }

    public double getRange() {
        return potentialAtomic.getRange();
    }

    public double energy(IMoleculeList molecules) {
        atoms.clear();
        for (int i=0; i<nBody; i++) {
            atoms.add(molecules.getMolecule(i).getChildList().getAtom(0));
        }
        return potentialAtomic.energy(atoms);
    }

    public void setBox(Box box) {
        potentialAtomic.setBox(box);
    }
}
