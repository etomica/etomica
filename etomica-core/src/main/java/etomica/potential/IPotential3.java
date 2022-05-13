/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.exception.MethodNotImplementedException;
import etomica.space.Vector;

/**
 * Methods for properties obtained for a soft, differentiable pair potential.
 *
 * @author David Kofke
 */
public interface IPotential3 {

    /**
     * Returns the range over which the potential applies.  IAtoms with a
     * greater separation do not interact.
     */
    default double getRange() {return Double.POSITIVE_INFINITY;}

    /**
     * The pair energy u(r^2) with no truncation applied.
     *  with r2xy the square of the distance between the particle x and y.
     */
    default double u(double r212, double r213, double r223) {
        throw new MethodNotImplementedException();
    }

    /**
     * Returns the energy between IAtoms atom1 and atom2 separated by vector
     * dr12 (which goes from atom1 to atom2; dr12 = r2 - r1).  PBC should not
     * be applied to dr12; PBC has already been accounted for and dr12 may even
     * exceed the Box dimensions when a lattice sum is used.  Likewise, the
     * IAtoms' actual positions (from getPosition()) ought to be ignored.
     */
    default double u(Vector dr12, Vector dr13, Vector dr23, IAtom atom1, IAtom atom2, IAtom atom3, double[] virial) {
        return u(dr12.squared(), dr13.squared(), dr23.squared());
    }

    default double udu(Vector dr12, Vector dr13, Vector dr23, IAtom atom1, IAtom atom2, IAtom atom3, double[] virial, Vector f1, Vector f2, Vector f3) {
        throw new MethodNotImplementedException();
    }
}