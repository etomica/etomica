/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Intermolecular potential that wraps an atomic pair potential.  The potential is
 * assumed to be monatomic.
 * 
 * @author Andrew Schultz
 */
public class P2MolecularMonatomic implements IPotentialMolecular {

    protected final Vector dr;
    protected final IPotential2 potentialAtomic;

    public P2MolecularMonatomic(Space space, IPotential2 potentialAtomic) {
        this.dr = space.makeVector();
        this.potentialAtomic = potentialAtomic;
    }

    public double energy(IMoleculeList molecules) {
        IAtom atom1 = molecules.get(0).getChildList().get(0);
        IAtom atom2 = molecules.get(1).getChildList().get(1);
        dr.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
        return potentialAtomic.u(dr, atom1, atom2);
    }
}
