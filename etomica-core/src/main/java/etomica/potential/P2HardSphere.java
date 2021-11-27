/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

/**
 * Basic hard-(rod/disk/sphere) potential.
 * Energy is infinite if spheres overlap, and is zero otherwise.  Collision diameter describes
 * size of spheres.
 * Suitable for use in space of any dimension.
 *
 * @author David Kofke
 */
public class P2HardSphere {

    public static P2HardGeneric makePotential(double sigma) {
        return new P2HardGeneric(new double[]{sigma}, new double[]{Double.POSITIVE_INFINITY});
    }
}//end of P2HardSphere
