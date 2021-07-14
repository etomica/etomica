/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.virial.mcmove;

import etomica.space.Space;
import etomica.util.random.IRandom;

public class MCMoveClusterAtomSQWTree extends MCMoveClusterAtomHSTree {

    protected final double pCore, pWell;

    public MCMoveClusterAtomSQWTree(IRandom random, Space _space, double lambda, double temperature) {
        super(random, _space, lambda);
        if (temperature == 0) {
            pCore = 0;
            pWell = 0;
            return;
        }
        double vCore = 1;
        double vWell = _space.powerD(lambda) - vCore;
        double Y = Math.exp(1 / temperature) - 1;
        if (Y > 1) {
            pCore = 0;
            pWell = vWell * (Y - 1) / (vWell + vCore + vWell * (Y - 1));
            return;
        }
        pCore = (1 - Y) * vCore / ((1 - Y) * vCore + Y * (vWell + vCore));
        pWell = 0;
    }

    protected double getSigma(int i, int j) {
        if (pCore > 0 && random.nextDouble() < pCore) {
            return 1;
        }
        if (pWell > 0 && random.nextDouble() < pWell) {
            return -sigma;
        }
        return sigma;
    }
}
