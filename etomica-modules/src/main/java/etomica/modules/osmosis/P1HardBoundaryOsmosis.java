/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.osmosis;

import etomica.atom.IAtomKinetic;
import etomica.box.Box;
import etomica.potential.P1HardBoundary;
import etomica.space.Vector;

/**
 * Potential that adds a wall at the center of the box (in x) on top of the
 * walls from P1HardBoundary.
 */
public class P1HardBoundaryOsmosis extends P1HardBoundary {

    public P1HardBoundaryOsmosis(Box box) {
        super(box.getSpace(), false, box);
    }

    @Override
    public double collisionTime(IAtomKinetic atom, Vector r, Vector v, int state, double falseTime) {
        double t = super.collisionTime(atom, r, v, state, falseTime);
        double x = r.getX(0);
        double vx = v.getX(0);
        double t0 = x * vx < 0 ? (-(x - Math.signum(x) * collisionRadius) / vx) : Double.POSITIVE_INFINITY;
        return t0 < 0 ? t : Math.min(t, t0);
    }

    @Override
    public int bump(IAtomKinetic atom, int oldState, Vector r, double falseTime, Vector deltaP, double[] du) {
        Vector v = atom.getVelocity();
        Vector dimensions = boundary.getBoxSize();
        double delmin = Double.MAX_VALUE;
        int imin = 0;
        //figure out which component is colliding
        for (int i = r.getD() - 1; i >= 0; i--) {
            double rx = r.getX(i);
            double vx = v.getX(i);
            double dxHalf = 0.5 * dimensions.getX(i);
            double del = (vx > 0.0) ? Math.abs(dxHalf - rx - collisionRadius) : Math.abs(-dxHalf - rx + collisionRadius);
            if (del < delmin) {
                delmin = del;
                imin = i;
            }
        }
        // check collision with central wall
        double del0 = Math.abs(r.getX(0)) - collisionRadius;
        if (Math.abs(del0) < delmin) {
            imin = 0;
        }

        deltaP.E(0);
        deltaP.setX(imin, -2 * v.getX(imin));
        deltaP.TE(atom.getType().getMass());
        v.setX(imin, -v.getX(imin));
        // dv = 2*NewVelocity
        double newP = atom.getPosition().getX(imin) - falseTime * v.getX(imin) * 2.0;
        atom.getPosition().setX(imin, newP);
        return 0;
    }
}
