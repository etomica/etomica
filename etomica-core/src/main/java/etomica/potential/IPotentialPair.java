/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.atom.AtomPair;
import etomica.atom.IAtom;
import etomica.space.Tensor;
import etomica.space.Vector;

public interface IPotentialPair extends IPotentialAtomic {

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
     * Computes the forces (and assigns them to f1 and f2) for IAtoms atom1
     * and atom2 separated by vector dr12 (which goes from atom1 to atom2;
     * dr12 = r2 - r1) and returns the energy between IAtoms atom1 and atom2.
     * PBC should not be applied to dr12; PBC has already been accounted for
     * and dr12 may even exceed the Box dimensions when a lattice sum is used.
     * Likewise, the IAtoms' actual positions (from getPosition()) ought to be
     * ignored.
     */
    double udu(Vector dr12, IAtom atom1, IAtom atom2, Vector f1, Vector f2);

    /**
     * Computes the forces (and adds them to f1 and f2) and torques (and
     * adds them to t1 and t2) for IAtoms atom1 and atom2 separated by
     * vector dr12 (which goes from atom1 to atom2; dr12 = r2 - r1) and returns
     * the energy between IAtoms atom1 and atom2.  PBC should not be applied to
     * dr12; PBC has already been accounted for and dr12 may even exceed the
     * Box dimensions when a lattice sum is used.  Likewise, the IAtoms' actual
     * positions (from getPosition()) ought to be ignored.
     */
    double uduTorque(Vector dr12, IAtom atom1, IAtom atom2, Vector f1, Vector f2, Vector t1, Vector t2);

    default Hessian d2u(Vector dr12, IAtom atom1, IAtom atom2) { return null; }

    /**
     * Integral used to evaluate correction to truncation of potential.
     */
    default double integral(double rC) {
        return 0;
    }

    default void u01TruncationCorrection(double[] uCorrection, double[] duCorrection) {

    }

    final class Hessian {
        public final Tensor r1r2, r1o2, r2o1, o1o1, o1o2, o2o2;
        public Hessian(Tensor r1r2, Tensor r1o2, Tensor r2o1, Tensor o1o1, Tensor o1o2, Tensor o2o2) {
            this.r1r2 = r1r2;
            this.r1o2 = r1o2;
            this.r2o1 = r2o1;
            this.o1o1 = o1o1;
            this.o1o2 = o1o2;
            this.o2o2 = o2o2;
        }

        public void PE(Hessian h2) {
            r1r2.PE(h2.r1r2);
            r1o2.PE(h2.r1o2);
            r2o1.PE(h2.r2o1);
            o1o1.PE(h2.o1o1);
            o1o2.PE(h2.o1o2);
            o2o2.PE(h2.o2o2);
        }
    }
}
