/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.space.Space;

/**
 * Factory that makes truncated and switched potentials.
 */
public class TruncationFactorySwitch implements TruncationFactory {

    protected final Space space;
    protected final double rc;
    protected final double switchFac;

    public TruncationFactorySwitch(Space space, double rc) {
        this(space, rc, 0.95);
    }

    public TruncationFactorySwitch(Space space, double rc, double switchFac) {
        this.space = space;
        this.rc = rc;
        this.switchFac = switchFac;
    }

    @Override
    public P2SoftSphericalSumTruncatedSwitched make(Potential2Soft... p2) {
        P2SoftSphericalSumTruncatedSwitched pTrunc = new P2SoftSphericalSumTruncatedSwitched(space, rc, p2);
        pTrunc.setSwitchFac(switchFac);
        return pTrunc;
    }
}
