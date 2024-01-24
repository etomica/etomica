/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.virial.mcmove;

import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Position source that returns a position in a sphere.  When constructed with a Boundary, the
 * displacement corresponds to a nearest image; the probability distribution is truncated at
 * the boundary.
 */
public class RandomPositionSphere implements MCMoveClusterAtomHSChain.InsertionPositionSource {

    public final Vector v;
    public final IRandom random;
    public final Boundary boundary;

    public RandomPositionSphere(Space space, IRandom random) {
        this(space, random, null);
    }

    public RandomPositionSphere(Space space, IRandom random, Boundary boundary) {
        v = space.makeVector();
        this.random = random;
        this.boundary = boundary;
    }

    public Vector position(double scale) {
        while (true) {
            v.setRandomInSphere(random);
            v.TE(scale);
            if (boundary==null) break;
            Vector x = boundary.centralImage(v);
            if (x.isZero()) break;
        }

        return v;
    }
}
