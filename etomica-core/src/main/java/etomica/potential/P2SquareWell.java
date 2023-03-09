/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

/**
 * Basic square-well potential.
 * Energy is infinite if spheres overlap, is -epsilon if less than lambda*sigma and not overlapping,
 * and is zero otherwise.  Core diameter describes size of hard core; lambda is multiplier to get range of well.
 * Suitable for use in space of any dimension.
 * Can be used with negative value for epsilon to produce square-shoulder potential. 
 *
 * This implementation determines which atoms are interacting (in the well) by
 * checking their distance.  This is prone to numerical roundoff issues, and so
 * atoms are bumped off the well after capture and escape.  This can cause
 * issues related to 3-particle collisions.  P2SquareWellRobust can be used
 * instead to avoid such issues.
 */
public class P2SquareWell {

    public static P2HardGeneric makePotential(double sigma, double lambda, double epsilon) {
        return new P2HardGeneric(new double[]{sigma, sigma * lambda}, new double[]{Double.POSITIVE_INFINITY, -epsilon});
    }
}
  
