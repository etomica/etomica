/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.virial.mcmove;

import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Position source that returns a position in a sphere.  When constructed with a Boundary, the
 * displacement corresponds to a nearest image; the probability distribution is truncated at
 * the boundary.
 */
public class RandomPositionLattice implements MCMoveClusterAtomHSChain.InsertionPositionSource {

    public final Vector v;
    public final IRandom random;
    public final int n2max;
    public final double L;

    public RandomPositionLattice(Space space, IRandom random, double L) {
        this(space, random, L, 1);
    }

    public RandomPositionLattice(Space space, IRandom random, double L, int n2max) {
        v = space.makeVector();
        this.random = random;
        this.n2max = 1;
        this.L = L;
    }

    public Vector position(double scale) {
        int n2;
        do {
            n2 = 0;
            for (int i=0; i<v.getD() && n2 <= n2max; i++) {
                int j = random.nextInt(L==2?2:3) - 1;
                v.setX(i, j);
                n2 += Math.abs(j);
            }
        } while (n2 > n2max);

        return v;
    }
}
