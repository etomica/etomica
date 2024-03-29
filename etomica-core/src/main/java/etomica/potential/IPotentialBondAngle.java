/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

/**
 * Interface for a potential capable of return bond angle energies and derivatives.
 */
public interface IPotentialBondAngle {

    /**
     * Returns the energy for the given value of cos(theta)
     */
    double u(double costheta);

    /**
     * Returns the energy for the given value of cos(theta),
     * r^2 between atoms 1 and 2 and r^2 between atoms 2 and 3
     */
    default double u(double costheta, double r12Sq, double r23Sq) {
        return u(costheta);
    }

    /**
     * Returns the energy (u) and du/dcostheta (du) as outparameters for the
     * given value of cos(theta)
     */
    void udu(double costheta, double[] u, double[] du);
}
