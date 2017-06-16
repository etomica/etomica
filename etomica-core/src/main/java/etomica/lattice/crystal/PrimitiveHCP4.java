/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.lattice.crystal;

import etomica.space.Space;

/**
 * Orthorhombic primitive cell for an HCP system.  Primitive-vector angles are
 * 90 degrees and the vectors are of unequal length (1, sqrt(3), sqrt(8/3)).
 * The hexagonal planes are perpendicular to the z axis.
 */
public class PrimitiveHCP4 extends PrimitiveOrthorhombic {

    /**
     * @param space better be a Space3D.  This is an IQ test.
     */
    public PrimitiveHCP4(Space space) {
        this(space, 1);
    }

    /**
     * @param space better be a Space3D.  This is an IQ test.
     * @param a the lattice constant -- vectors will have length
     *         a, a*sqrt(3), a*sqrt(8/3)
     */
    public PrimitiveHCP4(Space space, double a) {
        super(space, a, Math.sqrt(3)*a, Math.sqrt(8.0/3.0)*a);
    }
}
