/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

import etomica.space.Space;
import etomica.util.random.IRandom;

/**
 * Class that samples molecule positions based on a tree of hard spheres of
 * various sizes.  The interaction distance of each pair of spheres must be
 * given, which allows the composition and mixing rules to be set.
 * <p>
 * The work of handling trees is handled by the superclass.
 *
 * @author Andrew
 */
public class MCMoveClusterAtomHSTreeMix extends MCMoveClusterAtomHSTree {

    protected double[][] pairSigma;

    public MCMoveClusterAtomHSTreeMix(IRandom random, Space _space, double[][] pairSigma) {
        super(random, _space, 1);
        this.pairSigma = pairSigma;
    }

    public double getSigma(int i, int j) {
        return pairSigma[i][j];
    }
}
