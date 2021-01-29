/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.space.Space;

/**
 * Factory that makes truncated and force-shifted potentials.
 */
public class TruncationFactoryForceShift implements TruncationFactory {

    protected final Space space;
    protected final double rc;

    public TruncationFactoryForceShift(Space space, double rc) {
        this.space = space;
        this.rc = rc;
    }

    @Override
    public Potential2Soft make(Potential2SoftSpherical... p2) {
        return new P2SoftSphericalSumTruncatedForceShifted(space, rc, p2);
    }
}
