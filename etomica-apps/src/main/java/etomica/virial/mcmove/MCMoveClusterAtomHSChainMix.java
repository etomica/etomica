/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

import etomica.space.Space;
import etomica.util.random.IRandom;

/**
 * Class that samples molecule positions based on a chain of hard spheres of
 * various sizes.  The interaction distance of each pair of spheres must be
 * given, which allows the composition and mixing rules to be set.  The
 * sequence of spheres is chosen randomly.
 *
 * @author Andrew
 */
public class MCMoveClusterAtomHSChainMix extends MCMoveClusterAtomHSChain {

    protected double[][] pairSigma;

    public MCMoveClusterAtomHSChainMix(IRandom random, Space _space, double[][] pairSigma) {
        super(random, _space, 1);
        this.pairSigma = pairSigma;
    }

    public double getSigma(int i, int j) {
        return pairSigma[i][j];
    }
}
