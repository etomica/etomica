/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.sam;

import etomica.space.Vector;
import etomica.math.geometry.Plane;
import etomica.space.Space;

/**
 * Wrap a P1WCAWall and make it look like a Plane.  A boatload of
 * plane methods aren't overriden (there are a lot of them!) and calling them
 * will return garbage (or perhaps even crash).
 * DisplayBoxCanvasG3DSys only calls distanceTo and getD.
 *
 * @author Andrew Schultz
 */
public class WallPlane extends Plane {
    public WallPlane(Space space, P1WCAWall wallPotential) {
        super(space);
        this.wallPotential = wallPotential;
    }
    
    // DisplayBoxCanvasG3DSys calls this
    public double distanceTo(Vector v) {
        return v.getX(1) - wallPotential.getWallPosition()+2;
    }
    
    public double getA() {
        return 0;
    }
    
    public double getB() {
        return 1;
    }
    
    public double getC() {
        return 0;
    }

    // DisplayBoxCanvasG3DSys calls this
    public double getD() {
        return -wallPotential.getWallPosition()+2;
    }
    
    private static final long serialVersionUID = 1L;
    protected final P1WCAWall wallPotential;
}
