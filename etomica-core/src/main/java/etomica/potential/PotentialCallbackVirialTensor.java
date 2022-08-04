/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.potential.compute.PotentialCallback;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

public class PotentialCallbackVirialTensor implements PotentialCallback {

    protected final Space space;
    protected final Tensor virialTensor;

    public PotentialCallbackVirialTensor(Space space) {
        this.space = space;
        virialTensor = space.makeTensor();
    }

    public void reset() {
        virialTensor.E(0);
    }

    public Tensor getVirialTensor() {
        return virialTensor;
    }

    @Override
    public void pairCompute(int i, int j, Vector dr, double[] u012) {
        Vector f1 = space.makeVector();
        f1.Ea1Tv1(-u012[1] / dr.squared(), dr);
        virialTensor.PEv1v2(f1, dr);
    }
}
