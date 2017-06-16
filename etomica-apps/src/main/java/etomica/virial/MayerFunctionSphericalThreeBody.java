/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotentialAtomicMultibody;

/**
 * Non-additive Mayer function class for "spherical" potentials 
 *
 * @author Andrew Schultz
 */
public class MayerFunctionSphericalThreeBody extends MayerFunctionThreeBody {

    protected final IPotentialAtomicMultibody p3;
    
    public MayerFunctionSphericalThreeBody(IPotentialAtomicMultibody p3) {
        this.p3 = p3;
    }
    
    protected double energy(IMoleculeList molecules, double[] r2) {
        return p3.energy(r2);
    }

    public void setBox(Box box) {
        p3.setBox(box);
        super.setBox(box);
    }
}
