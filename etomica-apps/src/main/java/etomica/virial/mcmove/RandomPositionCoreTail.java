/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.virial.mcmove;

import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Position source that returns a position from a probability distrubtion consisting of a core
 * and tail.  When constructed with a Boundary, the displacement corresponds to a nearest image;
 * the probability distribution is truncated at the boundary.
 */
public class RandomPositionCoreTail implements MCMoveClusterAtomHSChain.InsertionPositionSource {

    public final Vector v;
    public final IRandom random;
    public final Boundary boundary;
    public final double pow, totPCore;
    public final int D;

    public RandomPositionCoreTail(Space space, IRandom random, double pow) {
        this(space, random, pow, null);
    }

    public RandomPositionCoreTail(Space space, IRandom random, double pow, Boundary boundary) {
        v = space.makeVector();
        this.random = random;
        this.boundary = boundary;
        this.pow = pow;
        D = space.D();
        totPCore = (pow-D) / pow;
    }

    public Vector position(double scale) {
        while (true) {
            double y = random.nextDouble();
            if (y < totPCore) {
                v.setRandomInSphere(random);
                v.TE(scale);
            }
            else {
                double r = Math.pow(1- (y - totPCore)/(1-totPCore), -1.0/(pow-D))*scale;
                v.setRandomSphere(random);
                v.TE(r);
            }
            if (boundary==null) break;
            Vector x = boundary.centralImage(v);
            if (x.isZero()) break;
        }

        return v;
    }
}
