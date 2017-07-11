/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.util.random.IRandom;
import etomica.space.Space;
import etomica.space.Vector;

public class RandomPositionSourceRectangular implements RandomPositionSource {

    protected final IRandom random;
    protected final Vector p;
    protected Box box;

    public RandomPositionSourceRectangular(Space space, IRandom random) {
        p = space.makeVector();
        this.random = random;
    }

    public Vector randomPosition() {
        p.setRandomCube(random);
        p.TE(box.getBoundary().getBoxSize());
        return p;
    }

    public void setBox(Box newBox) {
        box = newBox;
    }
}
