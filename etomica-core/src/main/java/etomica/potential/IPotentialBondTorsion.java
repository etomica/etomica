/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

/**
 * Interface for a potential capable of return bond torsion energies and derivatives.
 * Because it only takes cos(theta), the u(-theta) = u(theta) and u(0)=0.
 */
public interface IPotentialBondTorsion {

    /**
     * Returns the energy for the given value of cos(theta).
     */
    double u(double costheta);

    /**
     * Returns the energy (u) and du/d(cos(theta)) (du) as out parameters for the
     * given value of cos(theta).
     */
    void udu(double costheta, double[] u, double[] du);
}
