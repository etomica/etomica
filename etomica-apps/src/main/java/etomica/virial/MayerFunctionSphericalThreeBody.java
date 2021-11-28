/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.molecule.IMoleculeList;
import etomica.potential.Potential3Soft;

/**
 * Non-additive Mayer function class for "spherical" potentials 
 *
 * @author Andrew Schultz
 */
public class MayerFunctionSphericalThreeBody extends MayerFunctionThreeBody {

    protected final Potential3Soft p3;

    public MayerFunctionSphericalThreeBody(Potential3Soft p3) {
        this.p3 = p3;
    }

    protected double energy(IMoleculeList molecules, double[] r2) {
        return p3.u(r2[0], r2[1], r2[2]);
    }

    protected double energy(IMoleculeList molecules, double rAB2, double rAC2, double rBC2) {
        return p3.u(rAB2, rAC2, rBC2);
    }
}
