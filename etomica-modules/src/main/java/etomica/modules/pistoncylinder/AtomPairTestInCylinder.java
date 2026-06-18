/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.pistoncylinder;

import etomica.atom.AtomPairTest;
import etomica.atom.IAtom;
import etomica.space.Boundary;
import etomica.space.Vector;

public class AtomPairTestInCylinder implements AtomPairTest {

    public AtomPairTestInCylinder(Boundary boundary, P1HardMovingBoundary pistonPotential, double padding) {
        dimensions = boundary.getBoxSize();
        this.pistonPotential = pistonPotential;
        // bit flipper goes back and forth between 1 and 2
        bitFlipper = 1;
        this.padding = padding;
    }

    public int check(IAtom a) {
        double radius = pistonPotential.getCollisionRadius() + padding;

        Vector pos = a.getPosition();
        if (pos.getX(0) < -0.5 * dimensions.getX(0) + radius ||
                pos.getX(0) > 0.5 * dimensions.getX(0) - radius) {
            return 1;
        }
        if ((pos.getD() == 2 && (pos.getX(1) < pistonPotential.getWallPosition() + radius ||
                pos.getX(1) > 0.5 * dimensions.getX(1) - radius)) ||
                (pos.getD() == 3 && (pos.getX(1) > pistonPotential.getWallPosition() - radius ||
                        pos.getX(1) < -0.5 * dimensions.getX(1) + radius ||
                        pos.getX(2) < -0.5 * dimensions.getX(2) + radius ||
                        pos.getX(2) > 0.5 * dimensions.getX(2) - radius))) {
            return 1;
        }
        return 0;
    }


    public boolean test(IAtom atom1, IAtom atom2) {
        // always reject if both atoms are near a wall.  always accept if
        // both atoms are away from the wall.  If one is near and one not,
        // accept the pair half the time.  RDF needs this to avoid
        // over-counting pairs with one near the wall.  Ideally, we'd
        // accept them all and weight them half as much.
        int numOut = check(atom1) + check(atom2);
        // twiddle the last two bits, 1=>2, 2=>1
        // numOut=0 is always accepted, numOut=2 is never accepted
        bitFlipper ^= 3;
        return numOut < bitFlipper;
    }

    private final Vector dimensions;
    private final P1HardMovingBoundary pistonPotential;
    private int bitFlipper;
    private final double padding;
}
