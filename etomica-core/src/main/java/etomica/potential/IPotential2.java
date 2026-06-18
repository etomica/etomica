/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.exception.MethodNotImplementedException;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 * Methods for properties obtained for a soft, differentiable pair potential.
 *
 * @author David Kofke
 */
public interface IPotential2 {

    /**
     * Returns the range over which the potential applies.  IAtoms with a
     * greater separation do not interact.
     */
    default double getRange() {return Double.POSITIVE_INFINITY;}

    /**
     * The pair energy u(r^2) with no truncation applied.
     * @param r2 the square of the distance between the particles.
     */
    default double u(double r2) {
        throw new MethodNotImplementedException();
    }

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
    default double integral(Space space, double rC) {
        return 0;
    }

    default void u01TruncationCorrection(Space space, double[] uCorrection, double[] duCorrection) {

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

    /**
     * Returns the state of the pair of atoms (atom1 and atom2) at distance r12
     */
    default int getState(IAtom atom1, IAtom atom2, Vector r12) {
        throw new RuntimeException("Atom pairs have a state only for hard potentials");
    }

    /**
     * Returns the pair energy for the given state
     */
    default double getEnergyForState(int state) {
        throw new RuntimeException("Atom pairs have a state only for hard potentials");
    }

    /**
     * Computes the collision time for a given atom pair.
     *
     * @param atom1     the first atom
     * @param atom2     the second atom
     * @param r12       the position vector between the atoms (r2 - r1)
     * @param v12       the velocity difference between the atoms (v2 - v1)
     * @param state     the current state of the pair (overlapped, in well, etc).
     *                  the numbering scheme of states is an implementation
     *                  detail of the potential class, returned by getState or bump.
     * @param falseTime false-positioning time of the integration
     * @return the next collision time for a pair of atoms.
     */
    default double collisionTime(IAtomKinetic atom1, IAtomKinetic atom2, Vector r12, Vector v12, int state, double falseTime) {
        throw new RuntimeException("Atom pairs have a collision time only for hard potentials");
    }

    /**
     * Handles the collision between a pair of atoms.  The atom positions and velocities
     * are updated and the new state is returned.
     *
     * @param atom1     the first atom
     * @param atom2     the second atom
     * @param oldState  the state of the pair before the collision
     * @param r12       the vector between the pair at the collision (r2 - r1)
     * @param v12       the velocity difference between the atoms (v2 - v1)
     * @param falseTime the time into the future when the particles collide
     *                  the particle positions will be updated such that
     *                  r(falsetime) = r(0) + vnew * falsetime
     * @param virial    out parameter for the collision virial
     * @param du        the change in energy from this collision
     * @return the new state of the pair
     */
    default int bump(IAtomKinetic atom1, IAtomKinetic atom2, int oldState, Vector r12, Vector v12, double falseTime, double[] virial, double[] du) {
        throw new RuntimeException("Atom pairs can bump only for hard potentials");
    }

}