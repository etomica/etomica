/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

/**
 * Factory that makes (simple) truncated potentials.
 */
public class TruncationFactorySimple implements TruncationFactory {

    protected final double rc;

    public TruncationFactorySimple(double rc) {
        this.rc = rc;
    }

    @Override
    public IPotential2 make(IPotential2... p2) {
        return new P2SoftSphericalSumTruncated(rc, p2);
    }
}
