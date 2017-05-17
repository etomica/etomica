/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice;

import etomica.space.Space;

/**
 * Marker interface indicating that AbstractLattice.site(int[]) returns an
 * IVector instance. No specification is made here regarding whether repeated
 * calls will return the same IVector instance. However, repeated calls with the
 * same int[] argument should return vectors that are equal to each other
 * (assuming that the lattice is not modified between calls).
 */

public interface SpaceLattice extends AbstractLattice {

    public Space getSpace();

    public double[] getLatticeConstants();
}
