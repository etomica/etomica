/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.api.IBoundary;
import etomica.api.IRandom;
import etomica.api.IVector;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;

/**
 * Implementation of RandomPositionSource that can handle
 * BoundaryDeformablePeriodic.  If the boundary is not
 * BoundaryDeformablePeriodic, this class falls back to assuming rectangular
 * boundary.
 * 
 * @author Andrew Schultz
 */
public class RandomPositionSourceDeformable implements RandomPositionSource {
    
    public RandomPositionSourceDeformable(Space space, IRandom random) {
        p = space.makeVector();
        this.random = random;
    }

    public IVector randomPosition() {
        p.setRandomCube(random);
        IBoundary boundary = box.getBoundary();
        if (boundary instanceof BoundaryDeformablePeriodic) {
            ((BoundaryDeformablePeriodic)boundary).getBoundaryTensor().transform(p);
        }
        else {
            p.TE(boundary.getBoxSize());
        }
        return p;
    }

    public void setBox(Box newBox) {
        box = newBox;
    }

    protected final IRandom random;
    protected Box box;
    protected final IVector p;
}
