/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.exception.MethodNotImplementedException;
import etomica.space.Vector;

/**
 * Methods for properties obtained for a soft, differentiable pair potential.
 *
 * @author David Kofke
 */
public interface Potential2Soft extends Potential2Spherical, IPotentialPair {

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

    default Vector[][] gradientAndTorque(IAtomList atoms) {
        throw new MethodNotImplementedException();
    }

}