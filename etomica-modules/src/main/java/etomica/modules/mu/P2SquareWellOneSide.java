/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.mu;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.potential.P2HardGeneric;
import etomica.potential.Potential2Soft;
import etomica.space.Vector;

/**
 * Basic square-well potential.
 * Energy is infinite if spheres overlap, is -epsilon if less than lambda*sigma and not overlapping,
 * and is zero otherwise.  Core diameter describes size of hard core; lambda is multiplier to get range of well.
 * Suitable for use in space of any dimension.
 * Can be used with negative value for epsilon to produce square-shoulder potential.
 */
public class P2SquareWellOneSide extends P2HardGeneric implements Potential2Soft {

    public P2SquareWellOneSide(double coreDiameter, double lambda, double epsilon) {
        super(new double[]{coreDiameter, coreDiameter * lambda}, new double[]{Double.POSITIVE_INFINITY, -epsilon}, true);
    }

    @Override
    public int getState(IAtom atom1, IAtom atom2, Vector r12) {
        double x0 = atom1.getPosition().getX(0);
        double x1 = atom2.getPosition().getX(0);
        if (x0 < 0 || x1 < 0) {
            return -1;
        }
        return super.getState(atom1, atom2, r12);
    }

    @Override
    public double collisionTime(IAtomKinetic atom1, IAtomKinetic atom2, Vector r12, Vector v12, int collisionState, double falseTime) {
        double x0 = atom1.getPosition().getX(0) + atom1.getVelocity().getX(0) * falseTime;
        double x1 = atom2.getPosition().getX(0) + atom2.getVelocity().getX(0) * falseTime;
        if (x0 < 0 || x1 < 0) {
            return Double.POSITIVE_INFINITY;
        }
        return super.collisionTime(atom1, atom2, r12, v12, collisionState, falseTime);
    }

    @Override
    public double u(Vector dr12, IAtom atom1, IAtom atom2) {
        double x1 = atom1.getPosition().getX(0);
        double x2 = atom2.getPosition().getX(0);
        if (x1 < 0 || x2 < 0) {
            // on is ideal-gas
            return 0;
        }
        return super.u(dr12, atom1, atom2);
    }

    @Override
    public double energy(IAtomList pair) {
        throw new RuntimeException("nope");
    }

}
  
