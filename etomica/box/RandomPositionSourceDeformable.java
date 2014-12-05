/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.api.IVectorMutable;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;

/**
 * Implementation of RandomPositionSource that can handle
 * BoundaryDeformablePeriodic.  If the boundary is not
 * BoundaryDeformablePeriodic, this class falls back to assuming rectangular
 * boundary.
 * 
 * @author Andrew Schultz
 */
public class RandomPositionSourceDeformable implements RandomPositionSource {
    
    public RandomPositionSourceDeformable(ISpace space, IRandom random) {
        p = (IVectorRandom)space.makeVector();
        this.random = random;
    }

    public IVectorMutable randomPosition() {
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

    public void setBox(IBox newBox) {
        box = newBox;
    }

    protected final IRandom random;
    protected IBox box;
    protected final IVectorRandom p;
}
