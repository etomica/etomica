/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotential;
import etomica.potential.Potential2Spherical;

/**
 * @author kofke
 *
 * General Mayer function, which wraps the Mayer potential around an instance of
 * a Potential2 object.
 */
public class MayerGeneralSpherical implements MayerFunction {

    /**
     * Constructor Mayer function using given potential.
     */
    public MayerGeneralSpherical(Potential2Spherical potential) {
        this.potential = potential;
    }

    public double f(IMoleculeList pair, double r2, double beta) {
        double x = -beta*potential.u(r2);
        if (Math.abs(x) < 0.01) {
            return x + x*x/2.0 + x*x*x/6.0 + x*x*x*x/24.0 + x*x*x*x*x/120.0;
        }
        return Math.exp(x) - 1.0;
    }

    public IPotential getPotential() {
        return potential;
    }

    public void setBox(Box newBox) {
        potential.setBox(newBox);
    }

    private final Potential2Spherical potential;
}
