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

    default void udu(double r2, double[] u, double[] du) {
        u[0] = u(r2);
        du[0] = du(r2);
    }
}
