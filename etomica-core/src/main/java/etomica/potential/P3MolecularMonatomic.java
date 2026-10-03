/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Intermolecular 3-body potential that wraps an atomic 3-body  potential.  The potential is
 * assumed to be monatomic.
 * 
 * @author Andrew Schultz
 */
public class P3MolecularMonatomic implements IPotentialMolecular {

    protected final Vector dr12, dr13, dr23;
    protected final IPotential3 potentialAtomic;

    public P3MolecularMonatomic(Space space, IPotential3 potentialAtomic) {
        dr12 = space.makeVector();
        dr13 = space.makeVector();
        dr23 = space.makeVector();
        this.potentialAtomic = potentialAtomic;
    }

    public double energy(IMoleculeList molecules) {
        IAtom atom1 = molecules.get(0).getChildList().get(0);
        IAtom atom2 = molecules.get(1).getChildList().get(1);
        IAtom atom3 = molecules.get(1).getChildList().get(2);
        dr12.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
        dr13.Ev1Mv2(atom3.getPosition(), atom1.getPosition());
        dr23.Ev1Mv2(atom3.getPosition(), atom2.getPosition());
        return potentialAtomic.u(dr12, dr13, dr23, atom1, atom2, atom3, new double[1]);
    }
}
