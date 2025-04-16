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
public class RandomPositionLatticeSQW implements MCMoveClusterAtomHSChain.InsertionPositionSource {

    public final Vector v;
    public final IRandom random;
    public final int[][] degeneracy = new int[][]{{},{1,2},{1,4,4},{1,6,8,12}};
    public final int[][] degeneracyL2 = new int[][]{{},{1,1},{1,2,1},{1,3,3,1}};
    protected final double[] cprob;
    protected final int L;

    public RandomPositionLatticeSQW(Space space, IRandom random, double[] weights, int L) {
        v = space.makeVector();
        this.L = L;
        this.random = random;
        cprob = new double[weights.length];
        double sum = 0;
        int[] d = L==2 ? degeneracyL2[space.D()] : degeneracy[space.D()];
        for (int i=0; i<weights.length; i++) {
            sum += weights[i] * d[i];
        }
        cprob[0] = weights[0] * d[0] / sum;
        for (int i=1; i<weights.length; i++) {
            cprob[i] = cprob[i-1] + weights[i] * d[i] / sum;
        }
        cprob[cprob.length-1] = 1;
    }

    public Vector position(double scale) {
        double r = random.nextDouble();
        int in2;
        for (in2=0; in2<cprob.length; in2++) {
            if (cprob[in2] > r) break;
        }
        if (in2 == 0) {
            v.E(0);
            return v;
        }
        int n2;
        do {
            n2 = 0;
            for (int i=0; i<v.getD() && n2 <= in2; i++) {
                int j = random.nextInt(L==2 ? 2 : 3) - 1;
                v.setX(i, j);
                n2 += Math.abs(j);
            }
        } while (n2 != in2);

        return v;
    }
}
