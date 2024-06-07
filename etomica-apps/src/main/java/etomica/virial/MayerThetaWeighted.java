/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.potential.compute.PotentialCompute;

/**
 * MayerFunction implementation that returns f * du/dk where
 * u is intramolecular energy and k is some intramolecular parameter.
 * Optionally, U*exp(-beta*U) can be added, where U is the intermolecular energy.
 *
 * This class is useful when computing dbeta/dk at fixed B2 (B2=0 for theta point).
 */
public class MayerThetaWeighted implements MayerFunction {

    protected final MayerFunction f;
    protected PotentialCompute pcExtra;

    public MayerThetaWeighted(MayerFunction f) {
        this.f = f;
    }
    /**
     * Constructor Mayer function using given potential.
     */
    public void setPotentialExtra(PotentialCompute pcExtra) {
        this.pcExtra = pcExtra;
    }

    public double f(IMoleculeList pair, double r2, double beta) {
        return f.f(pair, r2, beta) * Math.exp(-beta*pcExtra.computeAll(false));
    }

    public void setBox(Box newBox) {
    }
}
