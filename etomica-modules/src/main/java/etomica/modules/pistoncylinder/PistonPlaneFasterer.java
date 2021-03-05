/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.pistoncylinder;

import etomica.math.geometry.Plane;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Wrap a P1HardMovingBoundary and make it look like a Plane.  A boatload of
 * plane methods aren't overriden (there are a lot of them!) and calling them
 * will return garbage (or perhaps even crash).
 * DisplayBoxCanvasG3DSys only calls distanceTo and getD.
 *
 * @author Andrew Schultz
 */
public class PistonPlaneFasterer extends Plane {
    public PistonPlaneFasterer(Space space, P1HardMovingBoundary pistonPotential) {
        super(space);
        this.pistonPotential = pistonPotential;
    }

    // DisplayBoxCanvasG3DSys calls this
    public double distanceTo(Vector v) {
        return v.getX(1) - pistonPotential.getWallPosition();
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
        return -pistonPotential.getWallPosition();
    }

    protected final P1HardMovingBoundary pistonPotential;
}
