/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.atom.AtomPair;
import etomica.atom.IAtom;
import etomica.space.Vector;

public interface IPotentialPair extends IPotentialTorque {

    /**
     * Returns the energy between IAtoms atom1 and atom2 separated by vector
     * dr12 (which goes from atom1 to atom2; dr12 = r2 - r1).  PBC should not
     * be applied to dr12; PBC has already been accounted for and dr12 may even
     * exceed the Box dimensions when a lattice sum is used.  Likewise, the
     * IAtoms' actual positions (from getPosition()) ought to be ignored.
     */
    default double u(Vector dr12, IAtom atom1, IAtom atom2) {
        return energy(new AtomPair(atom1, atom2));
    }

    /**
     * Computes the gradients (and assigns them to g1 and g2) for IAtoms atom1
     * and atom2 separated by vector dr12 (which goes from atom1 to atom2;
     * dr12 = r2 - r1) and returns the energy between IAtoms atom1 and atom2.
     * PBC should not be applied to dr12; PBC has already been accounted for
     * and dr12 may even exceed the Box dimensions when a lattice sum is used.
     * Likewise, the IAtoms' actual positions (from getPosition()) ought to be
     * ignored.
     */
    default double udu(Vector dr12, IAtom atom1, IAtom atom2, Vector g1, Vector g2) {
        return uduTorque(dr12, atom1, atom2, g1, g2, Vector.d(g1.getD()), Vector.d(g2.getD()));
    }

    /**
     * Computes the gradients (and assigns them to g1 and g2) and torques (and
     * assigns them to t1 and t2) for IAtoms atom1 and atom2 separated by
     * vector dr12 (which goes from atom1 to atom2; dr12 = r2 - r1) and returns
     * the energy between IAtoms atom1 and atom2.  PBC should not be applied to
     * dr12; PBC has already been accounted for and dr12 may even exceed the
     * Box dimensions when a lattice sum is used.  Likewise, the IAtoms' actual
     * positions (from getPosition()) ought to be ignored.
     */
    default double uduTorque(Vector dr12, IAtom atom1, IAtom atom2, Vector g1, Vector g2, Vector t1, Vector t2) {
        AtomPair pair = new AtomPair(atom1, atom2);
        double u = energy(pair);
        Vector[][] gt = gradientAndTorque(pair);
        g1.E(gt[0][0]);
        g2.E(gt[0][1]);
        t1.E(gt[1][0]);
        t2.E(gt[1][1]);
        return u;
    }
}
