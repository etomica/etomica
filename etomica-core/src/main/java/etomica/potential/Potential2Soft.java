/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;

/**
 * Methods for properties obtained for a soft, differentiable pair potential.
 *
 * @author David Kofke
 */
public interface Potential2Soft extends PotentialSoft, Potential2Spherical {

    /**
     * Hypervirial of the pair as given by the du(double) and d2u(double) methods
     */
    public double hyperVirial(IAtomList pair);

    /**
     * Integral used to evaluate correction to truncation of potential.
     */
    public double integral(double rC);

    /**
     * The derivative of the pair energy, times the separation r: r du/dr.
     */
    public double du(double r2);

    default void u012add(double r2, double[] u012) {
        u012[0] += u(r2);
        u012[1] += du(r2);
        u012[2] += d2u(r2);
    }

    /**
     * The second derivative of the pair energy, times the square of the
     * separation:  r^2 d^2u/dr^2.
     */
    double d2u(double r2);


    default void u01TruncationCorrection(double[] uCorrection, double[] duCorrection) {

    }
}
