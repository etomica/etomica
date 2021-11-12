/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.pistoncylinder;

import etomica.atom.IAtomList;
import etomica.atom.iterator.ApiLeafAtoms;
import etomica.atom.iterator.AtomsetIteratorBoxDependent;
import etomica.space.Boundary;
import etomica.space.Vector;


/**
 * Our own ApiFiltered that's box-dependent
 */
public class ApiFilteredCylinder extends ApiLeafAtoms implements AtomsetIteratorBoxDependent {
    public ApiFilteredCylinder(AtomFilterInCylinder filter) {
        super();
        this.filter = filter;
    }

    public IAtomList next() {
        IAtomList list = super.next();
        while (list != null && !filter.accept(list)) {
            list = super.next();
        }
        return list;
    }

    public int size() {
        int count = 0;
        reset();
        for (Object a = next(); a != null; a = next()) {
            count++;
        }
        return count;
    }


    protected final AtomFilterInCylinder filter;

    /**
     * Filter to expclude any pair with an atom within some distance from a
     * wall.
     */
    public static class AtomFilterInCylinder {
        public AtomFilterInCylinder(Boundary boundary, P1HardMovingBoundary pistonPotential, double padding) {
            dimensions = boundary.getBoxSize();
            this.pistonPotential = pistonPotential;
            this.padding = padding;
            // bit flipper goes back and forth between 1 and 2
            bitFlipper = 1;
        }

        public boolean accept(IAtomList atoms) {
            double radius = pistonPotential.getCollisionRadius() + padding;
            // always reject if both atoms are near a wall.  always accept if
            // both atoms are away from the wall.  If one is near and one not, 
            // accept the pair half the time.  RDF needs this to avoid 
            // over-counting pairs with one near the wall.  Ideally, we'd 
            // accept them all and weight them half as much. 
            int numOut = 0;
            for (int i = 0; i < 2; i++) {
                Vector pos = atoms.get(i).getPosition();

                if (pos.getX(0) < -0.5 * dimensions.getX(0) + radius ||
                        pos.getX(0) > 0.5 * dimensions.getX(0) - radius) {
                    numOut++;
                } else if ((pos.getD() == 2 && (pos.getX(1) < pistonPotential.getWallPosition() + radius ||
                        pos.getX(1) > 0.5 * dimensions.getX(1) - radius)) ||
                        (pos.getD() == 3 && (pos.getX(1) > pistonPotential.getWallPosition() - radius ||
                                pos.getX(1) < -0.5 * dimensions.getX(1) + radius ||
                                pos.getX(2) < -0.5 * dimensions.getX(2) + radius ||
                                pos.getX(2) > 0.5 * dimensions.getX(2) - radius))) {
                    numOut++;
                }
            }
            // twiddle the last two bits, 1=>2, 2=>1
            // numOut=0 is always accepted, numOut=2 is never accepted
            bitFlipper ^= 3;
            return numOut < bitFlipper;
        }

        private double padding;
        private final Vector dimensions;
        private final P1HardMovingBoundary pistonPotential;
        private int bitFlipper;
    }
}
