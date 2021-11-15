/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.exception.MethodNotImplementedException;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 * Methods for properties obtained for a soft, differentiable pair potential.
 *
 * @author David Kofke
 */
public interface Potential2Soft extends Potential2Spherical {

    /**
     * The derivative of the pair energy, times the separation r: r du/dr.
     */
    default double du(double r2) {
        throw new MethodNotImplementedException();
    }

    default void u012add(double r2, double[] u012) {
        u012[0] += u(r2);
        u012[1] += du(r2);
        u012[2] += d2u(r2);
    }

    /**
     * The second derivative of the pair energy, times the square of the
     * separation:  r^2 d^2u/dr^2.
     */
    default double d2u(double r2) {
        throw new MethodNotImplementedException();
    }

    /**
     * Returns the energy between IAtoms atom1 and atom2 separated by vector
     * dr12 (which goes from atom1 to atom2; dr12 = r2 - r1).  PBC should not
     * be applied to dr12; PBC has already been accounted for and dr12 may even
     * exceed the Box dimensions when a lattice sum is used.  Likewise, the
     * IAtoms' actual positions (from getPosition()) ought to be ignored.
     */
    default double u(Vector dr12, IAtom atom1, IAtom atom2) {
        return u(dr12.squared());
    }

    default double udu(Vector dr12, IAtom atom1, IAtom atom2, Vector f1, Vector f2) {
        double[] u012 = new double[3];
        double r2 = dr12.squared();
        u012add(r2, u012);
        Vector f = Vector.d(f1.getD());
        f.E(dr12);
        f.TE(u012[1] / r2);
        f1.PE(f);
        f2.ME(f);
        return u012[0];
    }

    default double uduTorque(Vector dr12, IAtom atom1, IAtom atom2, Vector f1, Vector f2, Vector t1, Vector t2) {
        return udu(dr12, atom1, atom2, f1, f2);
    }

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