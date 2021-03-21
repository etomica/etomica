/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

/**
 * A potential that returns the energy for a given electron density (used for
 * embedded atom method, EAM).
 */
public interface IPotentialEmbedding {

    /**
     * Returns the embedding energy for the given electron density rho.
     */
    double u(double rho);

    /**
     * Returns (as out parameters) the embedding energy (f) for the given
     * electron density rho along with the derivative of the energy with rho.
     */
    void udu(double rho, double[] f, double[] df);
}
